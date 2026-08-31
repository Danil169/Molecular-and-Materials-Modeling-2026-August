from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
import os

# Отключаем гипертрединг для чистого MPI
os.environ['OMP_NUM_THREADS'] = '1'

# 1. Структура кристалла кремния (Si)
atoms = Atoms(
    symbols=['Si']*2,
    positions=[
        [1.3574500000, 4.0723500000, 4.0723500000],
        [0.0000000000, 0.0000000000, 0.0000000000]
    ],
    cell=[
        [0.0000000000, 2.7149000000, 2.7149000000],
        [2.7149000000, 0.0000000000, 2.7149000000],
        [2.7149000000, 2.7149000000, 0.0000000000]
    ],
    pbc=[True, True, True]
)

# 2. Настройки запуска (Используем 4 ядра и правильный путь к псевдопотенциалам)
pw_command = 'mpirun -np 4 pw.x'
pw_profile = EspressoProfile(command=pw_command, pseudo_dir='../')
pseudopotentials = {'Si': 'Si.upf'}
scf_kpts = (2, 2, 2)

# 3. Список энергий обрезания (ecutwfc), которые мы хотим протестировать
cutoffs = [20, 30, 40, 50, 60, 70, 80]

print("Начинаем тест сходимости энергии обрезания (ecutwfc)...")
print(f"{'ecutwfc (Ry)':<15} | {'Полная Энергия (eV)':<20} | {'Разница (eV)':<15}")
print("-" * 55)

prev_energy = None

# Запускаем цикл расчетов
for ecut in cutoffs:
    input_data = {
        'control': {
            'calculation': 'scf',
            'pseudo_dir': '../'
        },
        'system': {
            'ibrav': 0, 'nat': 2, 'ntyp': 1,
            'ecutwfc': ecut,  # Меняем энергию обрезания в каждом цикле
            'occupations': 'smearing', 'smearing': 'gauss', 'degauss': 0.01
        },
        'electrons': {
            'conv_thr': 1.0e-10
        }
    }
    
    calc = Espresso(
        profile=pw_profile,
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=scf_kpts
    )
    
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    
    if prev_energy is None:
        de = 0.0
    else:
        de = energy - prev_energy
        
    print(f"{ecut:<15} | {energy:<20.6f} | {de:<15.6f}")
    prev_energy = energy

print("-" * 55)
print("Тест завершен!")
