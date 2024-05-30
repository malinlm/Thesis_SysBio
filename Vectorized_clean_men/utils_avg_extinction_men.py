#!/usr/bin/env python3
# Lanxin Zhang

from scipy import special
import numpy as np
from compute_pe_truvada import *


def p_inoculum(r, v0, vl=10001):
    """
    This function computes the probability of inoculum size v0 given a success rate r.
    parameters:
    r: float
        success rate of one exposure mode
    v0: int
        inoculum size of interest
    vl: int
        upper bound for donor viral load
    return:
    The probability of inoculum size v0 under the corresponding exposure mode. For details please check
    online Method.
    """
    mu = 4.51
    sigma = 0.98
    m = 0.3892

    def cdf(x):        # the CDF of donor viral load (log normal distribution)
        return (1 + special.erf((np.log10(x ** (1 / m)) - mu) / (2 ** 0.5 * sigma))) / 2

    f0 = np.arange(1, vl)
    import warnings
    warnings.filterwarnings("ignore")
    f_array = cdf(f0) - cdf(f0-1)
    row_mat, col_mat = np.meshgrid(v0, f0)
    matrix = special.comb(col_mat, row_mat) * (r ** row_mat) * ((1 - r) ** (col_mat - row_mat))
    return f_array @ matrix


def expected_pe(x, mode_dict):
    """
    Calculate the expected extinction probability for different exposure mode, according to the probability
    distribution of inoculum size.
    parameters:
    x: float
        the extinction probability of interest
    mode_dict: dict
        a dictionary which contains the probability (value) of corresponding inoculum size (key)
    """
    res = 0
    for k, v in mode_dict.items():
        res += v * x ** k
    return res


def process_pe_to_phi(indicator, pe_array): 
    """
    Process the results of compute_pe_truvada(indicator_arr), and generate the prophylactic efficacy.
    parameters:
    indicator: list
        boolean array with 4 elements, indicates the hypothesis (see function compute_pe_truvada)
    pe_array: list
        a list contains the results of function compute_pe_truvada (extinction probability of rai and rvi)

    """
    r_rvi = 9.1e-5  # success rate of rvi for inoculum size (see Supplementary text 2)
    r_rai = 3.7e-3  # success rate of rai
    dict_v2p_rvi = {i: p_inoculum(r_rvi, i)[0] for i in range(20)}     # probability (value) of each inoculum size (key)
    dict_v2p_rai = {i: p_inoculum(r_rai, i)[0] for i in range(20)}
    pe0 = 0.9017        # extinction probability without PrEP
    pe_rai_0 = expected_pe(pe0, dict_v2p_rai)
    pe_rvi_0 = expected_pe(pe0, dict_v2p_rvi)
    phi = None
    if indicator[3]:
        def mean_value_exposure_averaged(p_rai, p_rvi, ratio):  # calculate the averaged PE value based on RAI ratio
            return p_rai * ratio + (1 - ratio) * p_rvi

        ratio = 0.04    # ratio of rai/rvi among all sexual activities
        p0 = mean_value_exposure_averaged(pe_rai_0, pe_rvi_0, ratio)
        pe_rai = expected_pe(pe_array[0], dict_v2p_rai)
        pe_rvi = expected_pe(pe_array[1], dict_v2p_rvi)
        pe = mean_value_exposure_averaged(pe_rai, pe_rvi, ratio)
        phi = 1 - (1 - pe) / (1 - p0)
    elif indicator[4]: # lines 82-84 were added to account for 100% anal intercourse. 
        pe = expected_pe(pe_array[0], dict_v2p_rai)
        phi = 1 - (1 - pe) / (1 - pe_rai_0)
    else:
        pe = expected_pe(pe_array[0], dict_v2p_rvi)
        phi = 1 - (1 - pe) / (1 - pe_rvi_0)
    print(phi.shape)
    return phi