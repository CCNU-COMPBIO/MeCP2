# SASA and structural classification data

This directory contains the solvent-accessible surface area (SASA) data and derived structural classifications used to characterize MeCP2 MBD residues.

Structural analysis was performed using the MeCP2–DNA complex structure PDB 3C2I.

For each structurally resolved MBD residue, SASA values were calculated in the isolated-residue, protein-only, and protein–DNA complex contexts.

Residue exposure was characterized using:

Exposure ratio = SASA_protein / SASA_isolated

DNA burial was characterized using:

DNA-buried ratio = (SASA_protein - SASA_complex) / SASA_protein

Residues were classified according to the following criteria:

- DNA-buried ratio > 0.01: protein–DNA interface
- DNA-buried ratio ≤ 0.01 and exposure ratio < 0.20: protein interior
- DNA-buried ratio ≤ 0.01 and exposure ratio ≥ 0.20: protein surface

Protein–DNA interface classification was given priority over the exposure-based classification.

A solvent probe radius of 1.4 Å was used for SASA calculations.

These data were subsequently used to assign structurally resolved MeCP2 MBD missense variants to protein-interior, protein-surface, or protein–DNA interface categories.
