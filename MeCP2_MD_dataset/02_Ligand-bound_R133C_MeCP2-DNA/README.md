# Ligand-bound Arg133Cys MeCP2–DNA MD simulations

This directory contains ligand-bound molecular dynamics (MD) simulation data for Arg133Cys MeCP2–DNA complexes with five candidate compounds prioritized by virtual screening.

## Systems

- R133C_DB00266
- R133C_DB00445
- R133C_DB00694
- R133C_DB06237
- R133C_DB08794

## MD simulations

For each ligand-bound system, three independent 200-ns production MD simulations were performed using AMBER, corresponding to a cumulative simulation time of 600 ns per system.

The protein, DNA, and ligands were described using the ff14SB, OL15, and GAFF2 force fields, respectively, and the systems were solvated using the TIP4P-Ew water model.

The three independent production trajectories are denoted as:

- run1
- run2
- run3

## Analyses

The deposited data include:

- Protein–DNA complex backbone RMSD
- Protein Cα RMSF
- DNA-backbone RMSF
- Protein–DNA hydrogen bonds
- MM/GBSA estimates of MeCP2–DNA interaction energetics

RMSD was analyzed over 0–200 ns, whereas RMSF, hydrogen-bond, and MM/GBSA analyses used the 20–200 ns interval.

## Notes

DB00445, DB00694, and DB06237 remained associated with the Arg133Cys-centered MeCP2–DNA interface during the simulations. DB00266 and DB08794 did not maintain stable interface association and were not included in the detailed comparative analyses.
