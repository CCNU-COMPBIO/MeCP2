# Ligand-free MeCP2–DNA MD simulations

This directory contains ligand-free molecular dynamics (MD) simulation data for the WT and eight representative MeCP2 MBD mutant–DNA systems analyzed in this study.

## Systems

- WT
- R106G
- R106L
- R106W
- R111G
- R133C
- R133G
- R133L
- E143K

All systems were constructed based on the MeCP2 MBD–DNA complex derived from PDB entry 3C2I.

## MD simulations

For each system, three independent 500-ns production MD simulations were performed using AMBER, corresponding to a cumulative simulation time of 1.5 μs per system.

The protein and DNA were described using the ff14SB and OL15 force fields, respectively, and the systems were solvated using the TIP4P-Ew water model.

The three independent production trajectories are denoted as:

- run1
- run2
- run3

## Analyses

The deposited data include:

- Protein–DNA complex backbone RMSD
- Protein Cα RMSF
- DNA-backbone RMSF
- Intraprotein hydrogen bonds
- Protein–DNA hydrogen bonds
- MM/GBSA estimates of MeCP2–DNA interaction energetics

MeCP2 residues 94–162 were used consistently for the main trajectory analyses.

## Notes

These systems contain MeCP2 and DNA without bound small-molecule compounds. Ligand-bound Arg133Cys MeCP2–DNA simulations are provided separately in the `02_Ligand-bound_R133C_MeCP2-DNA` directory.
