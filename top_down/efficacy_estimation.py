import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from mpmath import power, mpf
from math import floor, ceil
from sampled_simulation import *
from collections import OrderedDict
import os
import seaborn as sns
import pandas as pd

def run_simulations_efficacy_uniform_sampled(n_tot_nodrug, n_tot_drug, py_nodrug, py_followup, n_inf_nodrug, n_simul=10):
    """
    sample efficacy from uniform distribution and run the simulation on group with detected drug
    parameters:
        n_tot_nodrug: total number of individuals of drug undetected subgroup
        n_tot_drug: total number of individuals of drug detected subgroup
        py_nodrug: total time of drug undetected subgroup  in person-years
        py_followup: total time in person-years
        n_inf_nodrug: number on infections in drug undetected subgroup
        n_siml: times of simulation performed
    Return: dictionary with number of infected as key; efficacy as value
    """
    print('Simulation for the drug-detected group initialized: efficacy sampled from uniform distribution')
    dict_infection2phi = dict()
    dict_infection2phi_forhistrogram = dict()
    t_total = 0
    phi = np.random.choice(np.arange(0, 1, 0.01), n_simul)
    pdf = compute_def_integral(n_tot_nodrug, py_nodrug, n_inf_nodrug)

    # infection rate sample from no drug group
    r_inf = np.random.choice([i/1e5 for i in range(len(pdf))], n_simul, p=pdf.astype(np.float64))
    #print("r_inf", r_inf)
    for i in range(n_simul):
             # simulate the sub-cohort of individuals with detectable drug.
             # simulate people who take the drug (with a sampled efficacy) and infection incidence taken from people with no drug (to see how efficacy changes this)
            t, total, infection = simulation_clinical_trial(n_tot_drug, py_followup, r_inf[i], phi[i])
            if infection[-1] not in dict_infection2phi:
                dict_infection2phi[infection[-1]] = []
            dict_infection2phi[infection[-1]].append(phi[i])
            dict_infection2phi_forhistrogram[infection[-1]] = phi[i]
            t_total += sum(t)
    # print('Average PY:', t_total/n_simul)  
    print('Simulation done')
    return dict_infection2phi, dict_infection2phi_forhistrogram 

def compute_efficacy_distribution(inf_dict, inf_nodrug_file, n_inf): 
    """
    process the data and compute the distribution of efficacy for each clinical trial
    parameters:
        inf_dict: dict containing the simulation data of infection numbers assuming a uniform distributed efficacy, 
                  generated using func: run_simulations_efficacy_uniform_sampled() 
        inf_nodrug_file: name of npy file containing the infection numbers of no-drug group in intervention arm, should be generated before
                         using stochastic simulation (here the data are pre-computed and stored in ../data/inf)
        n_inf: total number of infection in intervention arm for each clinical trial
    return:
        A dictionary contains the probability distribution of efficacy:  efficacy (key): frequency (value)
        Efficacy is discretely divided into 100 intervals, spanning from 0% to 100%.
    """
    folder_name = 'data'

    # Specify the file name
    file_name = inf_nodrug_file

    # Construct the full path to the file
    file_path = os.path.join(folder_name, file_name)

    data = np.load(file_path.format(inf_nodrug_file))
    #print("data")
    #print(data)
    p, infection = np.histogram(data, bins=np.unique(data), density=True)
    #print("p, infection")
    #print(p)
    #print(infection)

    #P(inf)
    infection_drug = n_inf - np.array(infection)

    #print("infection_drug")
    #print(infection_drug)
    n_inf_larger_than0 = sum(infection_drug < 0 )
    infection_drug = infection_drug[:-n_inf_larger_than0]
    p_drug = p[:-n_inf_larger_than0]
    p_drug= np.append(p_drug, 1-sum(p_drug))
    #print("p-drug")
    #print(p_drug)

    efficacy_dict = OrderedDict({i:0 for i in np.arange(0, 1, 0.01)}) # percentage efficacy between 0 and 100
    n_efficacy_total = sum([len(inf_dict[i]) for i in infection_drug if i in inf_dict])
    #print("n efficacy total")
    #print(n_efficacy_total)

    # for each infection number we have a probability distribution of efficacy
    print("for loop")
    #fig, ax = plt.subplots(floor(len(infection_drug)), ceil(len(infection_drug)), figsize=(10, 8))
    folder_name_figs = 'figures_P(effinf)_tdf2'

    for idx, inf in enumerate(infection_drug):
        #print(idx, inf)
        if inf in inf_dict:
            w1 = p_drug[idx]
            
            #this is the distribution for this infection
            
            phi, counts = np.unique(inf_dict[inf], return_counts=True)

            """
            plt.bar(phi*100, counts/len(inf_dict[inf]), label = f"infection number {inf}")
            plt.xlabel('Efficacy (%)')
            plt.ylabel('probability')
            plt.legend()
            file_path = os.path.join(folder_name_figs, f'P(eff_inf{inf}).png')
            plt.savefig(file_path)
            #plt.show()
            

            
            plt.plot(phi*100, counts/len(inf_dict[inf]), label = f"infection number {inf}")
            plt.xlabel('Efficacy (%)')
            plt.ylabel('probability')
            plt.legend()
            """
            for i, p in enumerate(phi):
                print(i,p)
                efficacy_dict[p] += (w1 * counts[i] / len(inf_dict[inf])) 
        else: print(inf)
    #plt.show()
    return efficacy_dict


if __name__ == "__main__":
        
    # IPREX
    """
    print("iprex")
    inf_dict_iprex, iprex_eff_sampled =run_simulations_efficacy_uniform_sampled(n_tot_nodrug=600, n_tot_drug=624, py_nodrug=720, py_followup=749, n_inf_nodrug=31, n_simul = 100000)
  
    # COMPUTE DISTRBUTION OF EFFICACY
    efficacy_dict_iprex = compute_efficacy_distribution(inf_dict = inf_dict_iprex, inf_nodrug_file = 'iPrEx.npy', n_inf = 34)

    # plot the efficacy vs probability
    x_values_iprex = []
    y_values_iprex = []
    for i in np.arange(0, 1, 0.01):
         x_values_iprex.append(i*100)
         y_values_iprex.append(efficacy_dict_iprex[i])


    dfiprex = pd.DataFrame(y_values_iprex)
    dfiprex.to_csv('efficacy_iprex.csv')
    

    ######################################################################################################
    """
    #HPTN083
    """
    print("hptn083")
    inf_dict_hptn083, hptn083_eff_sampled = run_simulations_efficacy_uniform_sampled(n_tot_nodrug=320, n_tot_drug=1964, py_nodrug=720, py_followup=2741, n_inf_nodrug=33, n_simul = 100000)

    efficacy_dict_hptn083 = compute_efficacy_distribution(inf_dict = inf_dict_hptn083, inf_nodrug_file = 'HPTN083_new3.npy', n_inf = 39)

    # plot the efficacy vs probability
    x_values_hptn083 = []
    y_values_hptn083 = []
    for i in np.arange(0, 1, 0.01):
         x_values_hptn083.append(i*100)
         y_values_hptn083.append(efficacy_dict_hptn083[i])

    dfhptn = pd.DataFrame(y_values_hptn083)
    dfhptn.to_csv('efficacy_hptn.csv')


    """
    
    
    # TDF2
    """
    print("tdf2")
    inf_dict_tdf2, td2_eff_sampled = run_simulations_efficacy_uniform_sampled(n_tot_nodrug=66, n_tot_drug=265, py_nodrug=85, py_followup=339, n_inf_nodrug=1, n_simul = 100000)

    efficacy_dict_tdf2 = compute_efficacy_distribution(inf_dict = inf_dict_tdf2, inf_nodrug_file = 'TDF2.npy', n_inf = 2)

    x_values_tdf2 = []
    y_values_tdf2 = []
    for i in np.arange(0, 1, 0.01):
         x_values_tdf2.append(i*100)
         y_values_tdf2.append(efficacy_dict_tdf2[i])

    dftdf2 = pd.DataFrame(y_values_tdf2)
    dftdf2.to_csv("efficacy_tdf2.csv")
    """
    # PARTNERS
    """
    print("partners")
    inf_dict_partners, partners_eff_sampled = run_simulations_efficacy_uniform_sampled(n_tot_nodrug=192, n_tot_drug=821, py_nodrug=317, py_followup=1355, n_inf_nodrug=3, n_simul = 100000)

    efficacy_dict_partners = compute_efficacy_distribution(inf_dict = inf_dict_partners, inf_nodrug_file = 'PARTNERS.npy', n_inf = 4)

    x_values_partners = []
    y_values_partners = []
    for i in np.arange(0, 1, 0.01):
         x_values_partners.append(i*100)
         y_values_partners.append(efficacy_dict_partners[i])

    dfpartners = pd.DataFrame(y_values_partners)
    dfpartners.to_csv("efficacy_partners.csv")
    """

    #IPERGAY
    """
    print("ipergay")
    inf_dict_ipergay, ipergay_eff_sampled = run_simulations_efficacy_uniform_sampled(n_tot_nodrug=18, n_tot_drug=171, py_nodrug=31, py_followup=189, n_inf_nodrug=2, n_simul = 100000)

    efficacy_dict_ipergay = compute_efficacy_distribution(inf_dict = inf_dict_ipergay, inf_nodrug_file = 'IPERGAY.npy', n_inf = 2)

    x_values_ipergay = []
    y_values_ipergay = []
    for i in np.arange(0, 1, 0.01):
         x_values_ipergay.append(i*100)
         y_values_ipergay.append(efficacy_dict_ipergay[i])

    dfipergay = pd.DataFrame(y_values_ipergay)
    dfipergay.to_csv("efficacy_ipergay.csv")
    """