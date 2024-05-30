# Model-based predictions of protective HIV PrEP levels in men who have sex with men
Master Thesis Systems Biology

There are different guidelines for PrEP (Pre-exposure Prophylaxis) for MSM (Men who have Sex with Men) and cisgender women. However, the clinical trials on which these guidelines are based have some limitations. Therefore, we aim to simulate these clinical trials to address these limitations and incorporate hypotheses aiming to explain adherence-efficacy relationships to determine if the efficacy is affected by any male/female differences. The code and data in this repository were used to answer the following research question and aims:

##### Research Question: Which concentration marker matches/correlates with the efficacy of PrEP in MSM?

To answer this research question, we will apply two independent approaches: a data-driven top-down approach and a mechanistic modelling bottom-up approach. These approaches are split into three aims:

Aim 1: To establish efficacy estimations for each of the trials using a data-driven, top-down approach via clinical trial simulations and compare incidence rates between heterosexual men and MSM.
Aim 2: To use mechanistic bottom-up modelling that incorporates factors hypothesized to influence efficacy (adherence, exposure site pharmacokinetics, dNTP levels, and exposure route) to simulate clinical trials and test the influence of concentration markers on prophylactic efficacy.
Aim 3: To use clinical top-down estimations to test the validity of each hypothesis in the bottom-up approach and identify the efficacy marker for MSM.

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
