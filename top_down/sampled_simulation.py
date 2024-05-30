import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from mpmath import power, mpf
from math import floor
import os

def compute_def_integral(n_tot, py, n_inf):
    """
    n_tot: initial population number
    n_inf: infection number
    py: total person-years
    return: a list of probability corresponding to infection incidence ranging in [0, n_tot/py] with step width 1e-5.
    """
    # the pdf
    def fun_to_integrate(x):
        return (power(x, n_inf)) * (power(n_tot - py * x , n_tot - n_inf))

    stepwidth = mpf(str(1e-5))  # use mpf for floats with higher precision

    #r_inf is in [0, n_tot/py]
    #small steps with which to integrate stepwise
    steps = floor(n_tot / py / stepwidth)
    #print(steps)

    #distribution will be given in the form of an array
    res = list()
    #x_ = list()
    #x = 0
    #integrate pdf of get cdf
    a_ = list()
    for i in range(steps):
        #x_.append(x)
        a = i * stepwidth
        b = (i + 1) * stepwidth
        fa = fun_to_integrate(a)
        fb = fun_to_integrate(b)
        res.append((fa + fb) * stepwidth / 2)
        a_.append(fa)

    res = np.array(res)
        #we want the normalized pdf
    # normalized
    res = res / res.sum()
    return res



def simulation_clinical_trial(n_individuals, py_followup, r_infection_incidence, phi=0):
    """
    Run simulation for clinical studies.
    For every individual enrolled, there's 1/avg_py probability to drop off and 1/infection_incidence
    probability to get infected.
    Parameters:
        n_individuals: total number of individuals
        py_followup: average followup years per person
        r_infection_incidence: infection rate
    Return: trajectories of time, total number of individuals and number of infections
    """

    np.random.seed()
    r_dropoff = n_individuals / py_followup - r_infection_incidence      # drop off rate in /person_year
    #print("1/py_followup", 1/py_followup )
    #print("r_dropoff", r_dropoff)
    #print("r_inf", r_infection_incidence)
    r_infection_incidence = r_infection_incidence * (1-phi) # 1-phi is for scaling infection incidence with sampled efficacy (phi)
    #print("r_inf", r_infection_incidence)
    t_list = [0] # list of timepoints
    n_total = [n_individuals] # number of individuals at each time point
    n_infection_list = [0]  # number of infections at each time point
    t = 0
    n_infection = 0
    while n_individuals > 0:
        # sum of reaction propensities to determine how likely the reactions occur per unit time
        # could also be seen as the total reaction rate
        # cumulative probability of all events
        # overall frequency of events in a system
        B = (r_dropoff + r_infection_incidence) * n_individuals
        # determine next reaction time by sampling from exponential distrbution
        tau = np.random.exponential(1/B)
        # update time
        t = t + tau
        # r is chosen r = U(0,1)
        r = np.random.random()
        # if reaction propensity of r_dropoff is greater than the threshold for making a stochastic decision (r*B)
        if r_dropoff * n_individuals > r * B:      # someone drop off
            n_individuals -= 1
        else:                      # someone got infected
            n_individuals -= 1
            n_infection += 1
        n_total.append(n_individuals)
        n_infection_list.append(n_infection)
        t_list.append(t)
    return t_list, n_total, n_infection_list


def clinical_simulation_incidence_sampled(n_individuals, py_followup, r_inf, n_simul=100000, phi=0):
    """
    Run simulation for one clinical trial with infection incidence sampled from the distribution computed
    by function 'compute_def_integral'. By default the simulation will be repeated 100000 times.
    Parameters:
        n_individuals: total number of individuals
        py_followup: average followup years per person
        r_inf: an array with length n_simul containing infection incidences sampled
    """
    t_mat = list()
    n_mat = list()
    inf_mat = list()
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals, py_followup, r_inf[i], phi)
        t_mat.append(t)
        n_mat.append(total)
        inf_mat.append(infection)
    return t_mat, inf_mat

if __name__ == "__main__":
    n_simul = 100000

    # PLACEBO rinf

    pdf_iprex = compute_def_integral(n_tot = 1217, py = 1460, n_inf = 64)
    r_inf_iprex = np.random.choice([i/1e5 for i in range(len(pdf_iprex))], n_simul, p=pdf_iprex.astype(np.float64))

    print('Drug undetected subgroup iPreX r_inf: ')
    print("Mean:", np.mean(r_inf_iprex), "median", np.median(r_inf_iprex), "95% quantile range: ", np.quantile(r_inf_iprex, (0.025, 0.975)))
    
    pdf_tdf2 = compute_def_integral(n_tot = 331, py = 424, n_inf = 10)
    r_inf_tdf2 = np.random.choice([i/1e5 for i in range(len(pdf_tdf2))], n_simul, p=pdf_tdf2.astype(np.float64))
    
    print('Drug undetected subgroup TDF2 r_inf: ')
    print("Mean:", np.mean(r_inf_tdf2), "median", np.median(r_inf_tdf2), "95% quantile range: ", np.quantile(r_inf_tdf2, (0.025, 0.975)))
    
    pdf_partners = compute_def_integral(n_tot = 963, py = 1589, n_inf = 24)
    r_inf_partners = np.random.choice([i/1e5 for i in range(len(pdf_partners))], n_simul, p=pdf_partners.astype(np.float64))

    print('Drug undetected subgroup PARTNERS r_inf: ')
    print("Mean:", np.mean(r_inf_partners), "median", np.median(r_inf_partners), "95% quantile range: ", np.quantile(r_inf_partners, (0.025, 0.975)))


    pdf_ipergay = compute_def_integral(n_tot = 201, py = 212, n_inf = 14)
    r_inf_ipergay = np.random.choice([i/1e5 for i in range(len(pdf_ipergay))], n_simul, p=pdf_ipergay.astype(np.float64))

    print('Drug undetected subgroup IPERGAY r_inf: ')
    print("Mean:", np.mean(r_inf_ipergay), "median", np.median(r_inf_ipergay), "95% quantile range: ", np.quantile(r_inf_ipergay, (0.025, 0.975)))

    # RUNNING SIMULATION ON DRUG UNDETECTED GROUPS

    # iPrEx
    """
    print("running iPrEx simulation")

    pdf_iprex = compute_def_integral(n_tot = 600, py = 720, n_inf = 31)
    r_inf_iprex = np.random.choice([i/1e5 for i in range(len(pdf_iprex))], n_simul, p=pdf_iprex.astype(np.float64))
    # r1 = u(0,1), r_inf = F^-1(r1)
    # [i/1e5 for i in range(len(pdf_iprex))] quantiles of our distribution (i.e. the x)
    # p=pdf_iprex.astype(np.float64) our CDF (between 0 and 1) (i.e. the u)
    # returns the smalles x s.t. F(x) >= u

    print('Drug undetected subgroup iPreX r_inf: ')
    print("Mean:", np.mean(r_inf_iprex), "median", np.median(r_inf_iprex), "95% quantile range: ", np.quantile(r_inf_iprex, (0.025, 0.975)))

    #print(r_inf_iprex)
    inf_list_iprex = []
    print("plotting simulation")
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals = 600, py_followup = 720, r_infection_incidence = r_inf_iprex[i], phi = 0)
        inf_list_iprex.append(infection[-1])
        plt.plot(t, infection)
    #t_mat, inf_mat = clinical_simulation_incidence_sampled(n_individuals = 600, py_followup = 720, r_inf = r_inf_iprex, n_simul=10, phi=0)
    np.array(inf_list_iprex)

    plt.xlabel('t')
    plt.ylabel('Number of Infected Individuals')
    plt.title('iPrEx infections in drug undetected subgroup')

    folder_name = "data"
    # Create the folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Save the array to a .npy file in the 'data' folder
    np.save(os.path.join(folder_name, 'iPrEx.npy'), inf_list_iprex)
    """
    
    
    # HPTN083
    """
    print("running HPTN083 simulation")

    pdf_hptn083 = compute_def_integral(n_tot = 320, py = 446, n_inf = 33)
    r_inf_hptn083 = np.random.choice([i/1e5 for i in range(len(pdf_hptn083))], n_simul, p=pdf_hptn083.astype(np.float64))
    
    print('Drug undetected subgroup HPNT 083 r_inf: ')
    print("Mean:", np.mean(r_inf_hptn083), "median", np.median(r_inf_hptn083), "95% quantile range: ", np.quantile(r_inf_hptn083, (0.025, 0.975)))
    
    #print(r_inf_iprex)
    inf_list_hptn083 = []
    print("infection list")
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals = 320, py_followup = 446, r_infection_incidence = r_inf_hptn083[i], phi = 0)
        inf_list_hptn083.append(infection[-1])
        # plt.plot(t, infection)
    
    np.array(inf_list_hptn083)

    folder_name = "data"
    # Create the folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Save the array to a .npy file in the 'data' folder
    np.save(os.path.join(folder_name, 'HPTN083_new6.npy'), inf_list_hptn083)
    
    print("file saved")
    """

    # TDF2
    """
    print("running TDF2 simulation")

    pdf_tdf2 = compute_def_integral(n_tot = 66, py = 85, n_inf = 1)
    r_inf_tdf2 = np.random.choice([i/1e5 for i in range(len(pdf_tdf2))], n_simul, p=pdf_tdf2.astype(np.float64))
    
    print('Drug undetected subgroup TDF2 r_inf: ')
    print("Mean:", np.mean(r_inf_tdf2), "median", np.median(r_inf_tdf2), "95% quantile range: ", np.quantile(r_inf_tdf2, (0.025, 0.975)))
    
    inf_list_tdf2 = []
    print("plotting simulation")
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals = 66, py_followup = 85, r_infection_incidence = r_inf_tdf2[i], phi = 0)
        inf_list_tdf2.append(infection[-1])

    folder_name = "data"

    np.save(os.path.join(folder_name, 'TDF2.npy'), inf_list_tdf2)
    
    print("file saved")
    """
    
    # PARTNERS
    """
    print("running PARTNERS simulation")

    pdf_partners = compute_def_integral(n_tot = 192, py = 317, n_inf = 3)
    r_inf_partners = np.random.choice([i/1e5 for i in range(len(pdf_partners))], n_simul, p=pdf_partners.astype(np.float64))

    print('Drug undetected subgroup PARTNERS r_inf: ')
    print("Mean:", np.mean(r_inf_partners), "median", np.median(r_inf_partners), "95% quantile range: ", np.quantile(r_inf_partners, (0.025, 0.975)))

    
    inf_list_partners = []
    print("plotting simulation")
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals = 192, py_followup = 317, r_infection_incidence = r_inf_partners[i], phi = 0)
        inf_list_partners.append(infection[-1])

    folder_name = "data"

    np.save(os.path.join(folder_name, 'PARTNERS.npy'), inf_list_partners)
    
    print("file saved")
    """

    # IPERGAY
    """
    print("running IPERGAY simulation")

    pdf_ipergay = compute_def_integral(n_tot = 18, py = 31, n_inf = 2)
    r_inf_ipergay = np.random.choice([i/1e5 for i in range(len(pdf_ipergay))], n_simul, p=pdf_ipergay.astype(np.float64))

    print('Drug undetected subgroup IPERGAY r_inf: ')
    print("Mean:", np.mean(r_inf_ipergay), "median", np.median(r_inf_ipergay), "95% quantile range: ", np.quantile(r_inf_ipergay, (0.025, 0.975)))

    
    inf_list_ipergay = []
    print("plotting simulation")
    for i in range(n_simul):
        t, total, infection = simulation_clinical_trial(n_individuals = 18, py_followup = 31, r_infection_incidence = r_inf_ipergay[i], phi = 0)
        inf_list_ipergay.append(infection[-1])

    folder_name = "data"

    np.save(os.path.join(folder_name, 'IPERGAY.npy'), inf_list_ipergay)
    
    print("file saved")
    """