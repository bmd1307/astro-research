import time
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.integrate as integrate

from galaxy_data import *
from supernova_data import *
from color_profile import *
from spitzer_profile import *
from rate_calc import *
from cross_reference import *
from tests import *


def sn_rate_by_radius(num_bins, sn_type, low_limit, high_limit):
    # parses the full sample of supernovae
    gal_dict = {}

    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

    # parses the color profiles
    color_dict = {}
    parse_color_profile_file(color_dict, 'nearby_galaxy_fuv_radius.txt')

    pair_galaxies_and_colors(gal_dict, color_dict)

    pair_galaxies_and_sne(gal_dict, sne_list)

    profiled_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.color_profile is not None]

    num_Ia = sum([curr_gal.sne_Ia for curr_gal in profiled_galaxies])
    num_SE = sum([curr_gal.sne_SE for curr_gal in profiled_galaxies])
    num_II = sum([curr_gal.sne_II for curr_gal in profiled_galaxies])

    print('SN counts for the fully profiled galaxies by type:')
    print('Ia:', num_Ia)
    print('SE:', num_SE)
    print('II:', num_II)

    # construct the CMFs for each of the galaxies
    for curr_gal in profiled_galaxies:
        curr_gal.construct_mass_profile()

    # calculates the limits for each of the bins with a low and high limit and a number of bins
    bin_limits = [low_limit + float(n) * (high_limit - low_limit) / num_bins for n in range(0, num_bins + 1)]

    sn_rate_bins = [0] * num_bins

    for bin_num in range(0, num_bins):
        curr_bin_low = bin_limits[bin_num]      # lower distance ratio limit for this bin
        curr_bin_high = bin_limits[bin_num + 1] # upper distance ratio limit for this bin
        curr_bin_sn_data = []

        for curr_gal in profiled_galaxies:
            #tuple structure: Num Sne, stellar mass, control time
            curr_range_num_sne = 0
            curr_range_mass = curr_gal.radial_range_mass(curr_bin_low, curr_bin_high)
            curr_range_control_time = None

            # count the number of SNe this galaxy has in the radial range
            for sn in curr_gal.supernovae:
                curr_sn_type = sn.sn_type_class()
                # convert the Ib, Ibc, Ic classifications to SE
                if curr_sn_type == 'Ib' or curr_sn_type == 'Ic' or curr_sn_type == 'Ibc':
                    curr_sn_type = 'SE'

                if curr_sn_type == sn_type and \
                        curr_bin_low < sn.distance_ratio('read') and sn.distance_ratio('read') < curr_bin_high:
                    curr_range_num_sne = curr_range_num_sne + 1

            if sn_type == 'Ia':
                curr_range_control_time = curr_gal.tc_Ia
            elif sn_type == 'SE':
                curr_range_control_time = curr_gal.tc_SE
            elif sn_type == 'II':
                curr_range_control_time = curr_gal.tc_II
            else:
                raise AttributeError(sn_type + ' is not a valid SN type: rate cannot be calculated')

            curr_bin_sn_data.append( (curr_range_num_sne, curr_range_mass, curr_range_control_time) )

        sn_rate_bins[bin_num] = sn_rate_general(curr_bin_sn_data)

        total_sne_for_this_bin = sum([data[0] for data in curr_bin_sn_data])
        print('Bin', bin_num, 'has', total_sne_for_this_bin, 'sne')

    print('Rates:')
    for bin_num in range(0, num_bins):
        print(bin_limits[bin_num], bin_limits[bin_num + 1], sn_rate_bins[bin_num] * 10**12) # convert to SNuM

    plt.plot(bin_limits[:-1], sn_rate_bins)
    plt.show()


def profile_params(mass_Msun):
    # returns sigma_0, h_r, sigma_bo, R_b, n

    log_stellar_mass = math.log10(mass_Msun)

    if 9.0 < log_stellar_mass < 9.5:
        return 51.1, 1.146, 178.6, 0.597, 1.681
    elif 9.5 < log_stellar_mass < 10.0:
        return 192.6, 1.975, 827, 0.171, 1.752
    elif 10.0 < log_stellar_mass < 10.5:
        return 500.3, 2.182, 4674.2, 0.167, 1.444
    elif 10.5 < log_stellar_mass < 11.0:
        return 733.9, 2.934, 24011.7, 0.101, 1.773
    else:
        return None


def mass_in_range_diaz_garcia(stellar_mass_Msun, r1_pc, r2_pc):
    params = profile_params(stellar_mass_Msun)
    if params is None:
        return None

    sigma_0, h_r, sigma_bo, R_b, n = params

    surface_density_func = lambda r_pc: sigma_bo * math.exp(-(r_pc * 0.001 / R_b)**(1/n)) + sigma_0 * math.exp(-r_pc * 0.001 / h_r)

    # integrate 2 pi r Sigma(r) from r1 to r2
    integral_result = integrate.quad((lambda r: 2 * math.pi * r * surface_density_func(r)), r1_pc, r2_pc)

    return integral_result[0]


def test_outskirts_mass_diaz_garcia():
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if 1e9 < (curr_gal.stellar_mass_Lum * 1e10) < 1e11]

    print('List of', len(list_galaxies), 'galaxies.')

    outskirts_ratio = lambda curr_gal: mass_in_range_diaz_garcia(curr_gal.stellar_mass_Lum * 1e10, curr_gal.r25_pc(), math.inf) /\
                                        mass_in_range_diaz_garcia(curr_gal.stellar_mass_Lum * 1e10, 0.0, math.inf)

    masses = [curr_gal.stellar_mass_Lum * 10**10 for curr_gal in list_galaxies]

    ratios = [outskirts_ratio(curr_gal) for curr_gal in list_galaxies]

    plt.scatter(masses, ratios)
    plt.xscale('log')
    plt.show()


def plot_num_sne_vs_radius():
    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    parse_sne_file(sne_list, 'sn-full-optimal.txt')

    g_centric_radii = [sn.distance_ratio("read") for sn in sne_list]

    g_centric_radii_1_plus = [curr_rad for curr_rad in g_centric_radii if curr_rad >= 1.0]

    g_centric_radii_1_to_2 = [curr_rad for curr_rad in g_centric_radii if curr_rad >= 1.0 and curr_rad <= 2.0]


    print('Total SNe                    ', len(g_centric_radii))
    print('SNe in outskirts (beyond r25)', len(g_centric_radii_1_plus))
    print('SNe from r25 < r < 2 * r25   ', len(g_centric_radii_1_to_2))

    print('50th percentile for outskirts SNe:', np.percentile(g_centric_radii_1_plus, 50))

    plt.hist(g_centric_radii, cumulative = True, normed = True, bins = freedman_diaconis_nbins(g_centric_radii))
    plt.show()

    plt.hist(g_centric_radii, bins = freedman_diaconis_nbins(g_centric_radii))

    plt.show()

    plt.hist(g_centric_radii, log = True, bins = freedman_diaconis_nbins(g_centric_radii))

    plt.show()

    plt.hist(g_centric_radii_1_plus, bins = freedman_diaconis_nbins(g_centric_radii_1_plus))

    plt.show()

    plt.hist(g_centric_radii_1_plus, bins = freedman_diaconis_nbins(g_centric_radii_1_plus), log = True)

    plt.show()

    plt.hist(g_centric_radii_1_to_2, log = True)

    plt.show()
