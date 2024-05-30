# Model-based predictions of protective HIV PrEP levels in men who have sex with men
Master Thesis Systems Biology

Code is informed by: https://github.com/KleistLab/PrEP_TruvadaWomen.git

This repository includes the code for recreating the top-down and bottom-up approaches.

#### top_down approach

sampled_simulation.py contains the functions for generating the sampled infection incidences and for running Gillespie's algorithm on each of the clinical trials. The simulation data can be found in folder "data" along with a function P_inf.py to calculate the probability distribution of infection. 

efficacy_estimation.py contains the functions for generating an efficacy distribution for each clinical trial by calling on sampled_simulation.py. The resulting efficacy distributions can be found in folder "efficacy_distribution".

hypothesis_tester.py utilizes the dataframes for probability of detectable plasma TFV "df_001.csv" and "df_035.csv" (two files containing the probability that TFV is detectable (above the lower limit of quantification, LLOQ) within different dosing adherence for LLOQ = 0.001uM and 0.035uM respectively)) in Vectorized_clean_men along with prophylactic efficacies in Vectorized_clean_men to simulate clinical trials and create new simulation data "HPTN083_vp.npy", "iPrEx_vp.npy", "IPERGAY_vp.npy".


#### Vectorized_clean_men

PEPV along with compute_pe_truvada.py and utils_avg_extinction men can be used to compute the PrEP efficacy trajectory of multiple regimens for multiple individuals in a single run. For a detailed example see example_men.ipynb.

bottom_up.ipynb is an notebook of the steps of the bottom_up approach. For simplicity and lower run time this is done using dosing.csv file and average PK parameters. The final p_value testing at the end requires the output of the P_inf.py as well as the output of hypothesis_tester.py in top_down.

MSM contains the pre-computed prophylactic efficacies for 1,000 virtual patients with different pk parameters (found in folder "pk_data").  The final p_value testing at the end requires the output of the P_inf.py as well as the output of hypothesis_tester.py in top_down.
