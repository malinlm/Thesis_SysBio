from sampled_simulation import *

import pandas as pd


def simulations_nodrug_and_drug_all_hypotheses(n_tot_nodrug, n_tot_drug, py_nodrug, py_followup, n_inf_nodrug, filename,
                                             df_conc, phi_dict, n_simul=10000):
    """
    Run simulation for different hypotheses (drug undetected subgroup and 8 hypotheses for drug detected group)
    and check the distribution of infection numbers. For each dosing scheme
    Parameters:
        n_tot_nodrug: total number of individuals of drug undetected subgroup
        n_tot_drug: total number of individuals of drug detected subgroup
        py_nodrug: total time of drug undetected subgroup  in person-years
        py_followup: total time in person-years
        n_inf_nodrug: number on infections in drug undetected subgroup
        filename: name of the file to save the results
        df_conc: a dataframe containing the probability of detectable plasma TFV (check supplementary text 1)
        phi_dict: dict containing the efficacy value of every hypothesis

    """
    pdf = compute_def_integral(n_tot_nodrug, py_nodrug, n_inf_nodrug)
    r_inf = np.random.choice([i / 1e5 for i in range(len(pdf))], n_simul, p=pdf.astype(np.float64))

    # drug undetected subgroup in intervention arm
    t, inf = clinical_simulation_incidence_sampled(n_tot_nodrug, py_nodrug, r_inf, n_simul)
    a = [inf[i][-1] for i in range(len(inf))]
    #print(a)
    py = [np.sum(t[i]) for i in range(len(t))]
    print('Drug undetected subgroup: ')
    print("average PYs: ", np.sum(py) / len(py), np.quantile(py, (0.025, 0.975)))
    print("Mean:", np.mean(a), "median", np.median(a), "95% quantile range: ", np.quantile(a, (0.025, 0.975)))
    print('-------------------------------------')
    res = [a]

    # drug detected subgroup in intervention arm

    weighted_phi = []
    print('Drug detected subgroup: ')
    for key, value in phi_dict.items():
        # calculate the weighted efficacy of 7 dosing schemes based on the percentage-concentration plots
        # formula for prophylactic efficacy
        # for iprex, tdf2, partners, ipergay
        # phi = sum([df_conc.iloc[i,0] * value[i] for i in range(7)]) / sum(df_conc.iloc[:,0])
        # for hptn 083
        phi = sum([df_conc[i] * value[i] for i in range(7)])
        weighted_phi.append(phi)
        print('Scenario indicator:', key, '  Average efficacy:', phi)
        t, inf = clinical_simulation_incidence_sampled(n_tot_drug, py_followup, r_inf, n_simul, phi=phi)
        py = [np.sum(t[i]) for i in range(len(t))]
        print("average PYs:", np.sum(py) / len(py), np.quantile(py, (0.025, 0.975)))
        a = [inf[i][-1] for i in range(len(inf))]
        print("Infections: mean:", np.mean(a), "95% quantile range:", np.quantile(a, (0.025, 0.975)))
        print('***')
        res.append(a)
    np.save('{}'.format(filename), np.array(res))
    print(weighted_phi)



if __name__ == "__main__":
    phi_dict = {'10001': [0.8256119498063447, 0.9062859807365935, 0.9563527894407503, 0.9661172635526987, 0.9713788263690433, 0.9716202007333552, 0.9717244186634107], 
                '10101': [0.8819727562353268, 0.9334627621631426, 0.9712051670220229, 0.971636617272747, 0.971860882345768, 0.971933411572823, 0.9719740111373857], 
                '11001': [0.7723918798266732, 0.8361740465470742, 0.8923314234662315, 0.8987491989417947, 0.9167563163292299, 0.9228908562111315, 0.9356077956347567], 
                '11101': [0.7903677472161225, 0.8474998986958393, 0.9017957828137765, 0.9069416100158988, 0.9243532704278015, 0.9300317097607239, 0.9423560577518774]}
    
    phi_dict2 = {'10001': [0.5743758, 0.7996121, 0.90017134, 0.9313083, 0.95989025, 0.97156006, 0.9784501,], 
                '10101': [0.6503719, 0.8614176, 0.94225436, 0.963688, 0.9816999, 0.98692375, 0.98902905], 
                '11001': [0.418226812, 0.580964383, 0.671698391, 0.719404577, 0.760564603, 0.789913984, 0.812941811], 
                '11101': [0.45383775, 0.6224243, 0.71343863, 0.7583164, 0.7973944, 0.82406145, 0.84467405]}
    

    # plot the phi_dict
    mean_vals = []
    key_vals = []
    for key in phi_dict.keys():
        key_vals.append(key)
        m = np.mean(np.asarray(phi_dict[key]))
        mean_vals.append(m)
        plt.scatter(range(1,8), np.asarray(phi_dict[key])*100, 
                    marker = "x", linewidth = 2, label = f"{key}")
        plt.plot(range(1,8), np.asarray(phi_dict[key])*100, linewidth = 1.2)

    plt.plot(range(1,8), [50]*7, linestyle = "--", color = "red", linewidth = 0.7)
    plt.plot(range(1,8), [90]*7, linestyle = "--", color = "darkred",linewidth = 0.7)
    plt.ylim((0,110))
    plt.xlabel("Dose per week")
    plt.ylabel("Prophylactic efficacy (%)")
    plt.legend()
    plt.show()
    

    print("")
    print("mean efficacy values for bottom-up hypotheses")
    print(key_vals)
    print(mean_vals)

    
    df_001 = pd.read_csv('Vectorized_clean_men/df_001.csv', index_col=0)  # LLOQ = 0.001uM IPERGAY
    df_035 = pd.read_csv('Vectorized_clean_men/df_035.csv', index_col=0)  # LLOQ = 0.035uM iPrEx

    df_001 = df_001/100
    df_035 = df_035/100
    # n_tot_nodrug, n_tot_drug, py_nodrug, py_followup, n_inf_nodrug, filename, df_conc, phi_dict, n_simul=1000
    df_hptn = [0.0956, 0.054, 0.054, 0.1411, 0.1411, 0.1411, 0.3676]
    
    #n_tot_nodrug=320, n_tot_drug=1964, py_nodrug=720, py_followup=2741, n_inf_nodrug=37
    print(df_001)
    #print(df_001[0].iloc[0])  
    print("") 
    print(df_001['0'])
    for key, value in phi_dict.items():
        # calculate the weighted efficacy of 7 dosing schemes based on the percentage-concentration plots
        phi = sum([df_001.iloc[i,0] * value[i] for i in range(7)]) / sum(df_001.iloc[:,0])
        print(phi)
    
    # HPTN 083
    print("HPTN 083")
    simulations_nodrug_and_drug_all_hypotheses(n_tot_nodrug=320, n_tot_drug=1964, py_nodrug=446, py_followup=2741, n_inf_nodrug=33,
                                                filename = "HPTN083_vp", df_conc = df_hptn, phi_dict = phi_dict2, n_simul = 10000)
    """ 
    # iPrEx
    print("iPrEx")
    simulations_nodrug_and_drug_all_hypotheses(n_tot_nodrug=600, n_tot_drug=624, py_nodrug=720, py_followup=749, n_inf_nodrug=31,
                                                filename = "iPrEx_vp", df_conc = df_035, phi_dict = phi_dict2, n_simul = 10000)
    #""" 
                                               
    """ 
    print("IPERGAY")
    simulations_nodrug_and_drug_all_hypotheses(n_tot_nodrug=18, n_tot_drug=171, py_nodrug=31, py_followup=189, n_inf_nodrug=2,
                                                filename = "IPERGAY_vp", df_conc = df_001, phi_dict = phi_dict2, n_simul = 10000)
    """
   