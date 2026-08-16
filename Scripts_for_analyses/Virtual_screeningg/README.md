# Virtual screening

This directory contains the scripts and configuration files used for the structure-based virtual-screening workflow targeting the Arg133Cys-centered MeCP2–DNA interface.

The workflow consisted of two stages:

1. Docking against the intact Arg133Cys MeCP2–DNA complex.
2. Counter-screening against a DNA-removed model generated from the same representative structure.

Ligands were prepared from the approved-compound library obtained from DrugBank.

AutoDock Vina was used for molecular docking.

Compounds with favorable docking scores at the intact MeCP2–DNA interface were further evaluated using the DNA-removed model. Compounds with a docking score ≤ -7 kcal/mol in the DNA-removed model were deprioritized.

Docking scores were used for relative prioritization and were not interpreted as experimental binding free energies.
