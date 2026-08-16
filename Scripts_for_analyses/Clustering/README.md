# Trajectory clustering

This directory contains the CPPTRAJ input and analysis scripts used for clustering the Arg133Cys MeCP2–DNA MD trajectories.

Frames from the three independent Arg133Cys production trajectories were combined for clustering.

The protein backbone of MeCP2 residues 94–162 together with the DNA backbone was used for structural alignment and RMSD-based clustering.

K-means clustering was evaluated for cluster numbers from K = 2 to K = 10.

Cluster quality was assessed using the Davies–Bouldin index (DBI) and pseudo-F statistic (pSF).

The K = 2 solution was selected, and the centroid of the dominant cluster was used as the representative structure for subsequent virtual screening.
