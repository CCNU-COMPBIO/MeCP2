# SASA and structural classification

This directory contains the scripts used for solvent-accessible surface area (SASA) calculations and structural classification of MeCP2 MBD residues using the MeCP2–DNA structure PDB 3C2I.

Residue exposure was characterized using:

Exposure ratio = SASA_protein / SASA_isolated

DNA burial was characterized using:

DNA-buried ratio = (SASA_protein - SASA_complex) / SASA_protein

Residues with a DNA-buried ratio > 0.01 were classified as protein–DNA interface residues.

For non-interface residues:

- Exposure ratio < 0.20: protein interior
- Exposure ratio ≥ 0.20: protein surface

A solvent probe radius of 1.4 Å was used for the SASA calculations.
