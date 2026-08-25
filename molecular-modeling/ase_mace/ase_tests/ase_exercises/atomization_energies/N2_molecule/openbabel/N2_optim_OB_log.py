import os
from openbabel import openbabel as ob
from openbabel import pybel

def run_benchmarks_with_ff_metadata(output_log="n2_optimization_summary.txt", geom_dir="n2_optimized_geometries"):
    # 1. Setup Environment
    if not os.path.exists(geom_dir):
        os.makedirs(geom_dir)

    # Initial N2: Stretched to 1.25 A
    d_initial = 1.25
    # The second line is the "Comment" line where we store the FF name
    xyz_template = "2\nN2 optimized with {ff_name}\nN 0.0 0.0 0.0\nN 0.0 0.0 {dist}"

    # Configure OB Logger
    ob_log = ob.obErrorLog
    ob_log.SetOutputLevel(4)

    print(f"Starting benchmarks. Saving files...")

    with open(output_log, "w") as log_file:
        log_file.write("N2 Optimization Benchmark\n")
        log_file.write("="*60 + "\n\n")

        for ff_name in pybel.forcefields:
            # Create fresh molecule using the template
            initial_xyz = xyz_template.format(ff_name=ff_name, dist=d_initial)
            mol = pybel.readstring("xyz", initial_xyz)
            ff = pybel._forcefields[ff_name]

            if not ff.Setup(mol.OBMol):
                log_file.write(f"FORCE FIELD: {ff_name.upper()} - SETUP FAILED\n\n")
                continue

            # 2. Optimize
            ob_log.ClearLog()
            ff.ConjugateGradients(500, 1.0e-6)
            ff.GetCoordinates(mol.OBMol)

            # 3. Data Extraction
            a1, a2 = mol.atoms[0].coords, mol.atoms[1].coords
            final_dist = ((a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2)**0.5
            energy = ff.Energy()

            # 4. Update the molecule's title so the exported XYZ has the FF name
            mol.title = f"N2_optimized_with_{ff_name}"
            
            geom_filename = os.path.join(geom_dir, f"n2_opt_{ff_name}.xyz")
            mol.write("xyz", geom_filename, overwrite=True)

            # 5. Log writing
            log_file.write(f"FORCE FIELD: {ff_name.upper()}\n")
            log_file.write(f"Final Bond Length: {final_dist:.4f} Angstroms\n")
            log_file.write(f"Final Energy:      {energy:.4f} kJ/mol\n")
            
            messages = ob_log.GetMessagesOfLevel(0)
            if messages:
                log_file.write("--- Internal Debug Logs ---\n")
                for msg in messages:
                    log_file.write(f"  [OB]: {msg.strip()}\n")
            
            log_file.write("-" * 60 + "\n\n")
            ob_log.ClearLog()

    print(f"Done! Check the '{geom_dir}' folder for the XYZ files.")

if __name__ == "__main__":
    run_benchmarks_with_ff_metadata()

