import numpy as np
from random import uniform
import matplotlib.pyplot as plt
import os

def plot_p_inf(inf_nodrug_file, n_inf):

    folder_name = 'data'

    # Specify the file name
    file_name = inf_nodrug_file

    # Construct the full path to the file
    file_path = os.path.join(folder_name, file_name)

    data = np.load(file_path.format(inf_nodrug_file))
    #print("data")
    #print(data)
    p, infection = np.histogram(data, bins=np.unique(data), density=True)
    print(p,infection)
    print(f'Drug undetected subgroup {inf_nodrug_file} infection numbers: ')
    print("Mean:", np.mean(data), "median", np.median(data), "95% quantile range: ", np.quantile(data, (0.025, 0.975)))

    #plt.plot(infection, p)
    plt.bar(np.unique(data)[:-1], p, width=np.diff(np.unique(data)), align='edge', alpha=0.5, label='p')
    #plt.xlabel('Number of Infected Individuals in drug undetected cohort')
    plt.ylabel('Probability', fontdict={'family': 'serif', 'size': 14})
    plt.ylim(0, 0.3)
    plt.xlim(0,65)
    plt.title(f'{inf_nodrug_file} P(inf)')
    plt.show()

    infection_drug = n_inf - np.array(infection)

    n_inf_larger_than0 = sum(infection_drug < 0 )
    infection_drug = infection_drug[:-n_inf_larger_than0]
    p_drug = p[:-n_inf_larger_than0]
    p_drug= np.append(p_drug, 1-sum(p_drug))
    p_inf_dict = dict(zip(infection_drug, p_drug))

    print("Drug detected group probabilites")
    print(p_inf_dict)
    plt.bar(infection_drug, p_drug)
    #plt.xlabel('Number of Infected Individuals in drug detected cohort')
    plt.ylabel('Probability', fontdict={'family': 'serif', 'size': 14})
    plt.title(f'{inf_nodrug_file} P(inf)')
    plt.ylim(0, 0.5)
    plt.xlim(0, 30)
    plt.show()

print("iprex")
plot_p_inf(inf_nodrug_file = 'iPrEx.npy', n_inf = 34)

print('hptn083')
plot_p_inf(inf_nodrug_file = 'HPTN083_new6.npy', n_inf = 39)

print('tdf2')
plot_p_inf(inf_nodrug_file = 'TDF2.npy', n_inf = 2)

print('partners')
plot_p_inf(inf_nodrug_file = 'PARTNERS.npy', n_inf = 4)

print('ipergay')
plot_p_inf(inf_nodrug_file = 'IPERGAY.npy', n_inf = 2)

