===========================================
Molecular and Materials Modeling — Final Package
===========================================

This repository contains the completed work for the Two-Day Hybrid School on
Molecular and Materials Modeling (BLTP JINR, Dubna; 26–27 August 2026).

Final reports
=============

The current, submission-ready reports are kept in the repository root:

* ``Report_Materials_Modeling_EN.docx`` — English report; submit this file to
  the instructor.
* ``Report_Materials_Modeling_RU.docx`` — Russian counterpart for reference.

Both reports contain Part I (molecular modeling) and Part II (materials
modeling), the eight numbered figures, and the corrected OpenXML table
structure. The reports were validated after the final revision.

Repository layout
=================

``molecular-modeling/``
    Part I: software checks, performance benchmarks, ASE/MACE, molecular
    properties, U-complex, VSEPR, ionization-potential, and GROMACS exercises.

``materials-modeling/``
    Part II: Quantum ESPRESSO calculations for Si, Al, graphene, work function,
    and adsorption. ``report_figures_part1/`` contains the two Part I figures
    embedded in the final reports.

``Teaching-materials-2026August.docx``
    Original course handout.

``create_reports.py`` and ``update_report_part1.py``
    Report-generation and final-update provenance scripts. They are retained as
    supporting workflow material, not as the submission files.

Using the repository
====================

Each Git-tracked directory contains a ``readme.rst`` describing its purpose and
the associated inputs, scripts, and results. Calculation artifacts are retained
for reproducibility; no calculation needs to be rerun to read the submitted
report.

The original course page is https://indico.jinr.ru/event/6303/.
