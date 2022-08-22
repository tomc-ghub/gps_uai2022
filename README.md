# gps_uai2022
README file for GPS Matlab code used in article ‘Greedy Equivalence Search in the Presence of Latent Confounders’ (Claassen & Bucur, UAI2022)

Main script:
=> Test_runGPS.m : example wrapper script for running GPS by calling Run_GPS.m
 Configure parameters:
 - GPS run type:  select GPS version baseline/extended/hybrid (see article)
 - ScoreType : scoring metric ‘MML’ (Gaussian MAG score with BIC penalty) or ‘SHD’ (structural Hamming distance)
 - choose example test case (graph_id=2-5) or random generated graph / MVG model (graph_id=1)
 - for details on other parameters see description/comments in file

Other key scripts:
- Run_GPS.m : core GPS routine that runs the greedy equivalence search for the chosen scoring metric (Algorithm 4 in article)
- GetNeighbourMECs.m : function that generates a collection of PAG/MEC neighbours for a given PAG (Definition 2, Algorithm 3)
- mag_to_mec.m : converts maximal/partial ancestral graph into MEC representation (Algorithm 1)
- mec_to_cpag.m : converts MEC into corresponding completed Partial Ancestral Graph (Algorithm 2)
- experimental data generated using 'mk_ag.m' (random ancestral graph) and 'ag_to_random_MVG_pdf.m' (matching multivariate Gaussian model)
