# Brennan Dell
# Supernova Rate in Outskirts

import math
import matplotlib.pyplot as plt
import numpy as np

import galaxy_data
import supernova_data


def pair_galaxies_and_sne(gal_dict, sne_list):
    print('Pairing each the galaxies from Graur with Leaman data')
    hostless_sne_count = 0

    for gal in gal_dict.values():
        gal.supernovae = []

    for curr_sn in sne_list:
        g_name = curr_sn.galaxy_name
        if g_name in gal_dict.keys():
            gal_dict[g_name].supernovae.append(curr_sn)
        else:
            hostless_sne_count = hostless_sne_count + 1

    if hostless_sne_count > 0:
        print(hostless_sne_count, 'SNe were not paired with their host galaxy')
    else:
        print('All SNe paired with their host galaxy.')


# normalizes each value in a list to 0.0 - 1.0
# ex [1, 2, 3, 5] -> [0, 0.25, 0.5, 1.0]
def normalize(arr):

    max_val = max(arr)
    min_val = min(arr)

    return [(n - min_val) / (max_val - min_val) for n in arr]


# Recreates the plot in Figure 8 of Leaman et al., the cumulative distributions of the SN types (Ia, Ib, Ibc, Ic, II)
def plot_sne_cumulative():

    # parses the full sample of supernovae
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full.txt')

    # remove the galaxies with a error in the position angle (-99.999)
    num_all_sne = len(sne_list)
    sne_list = [sn for sn in sne_list if sn.galaxy_pa_deg >= 0.0]
    print('Removed', str(num_all_sne - len(sne_list)), ' SNe with host galaxies that had bad position angles')

    # create a new list for each type of supernovae

    # Ia list
    sne_Ia = [sn for sn in sne_list if sn.sn_type_class() == 'Ia']
    print('Found', len(sne_Ia), 'SNe of type Ia')

    # Ib list
    sne_Ib = [sn for sn in sne_list if sn.sn_type_class() == 'Ib']
    print('Found', len(sne_Ib), 'SNe of type Ib')

    # Ibc list (is actually the combination of Ib, Ibc and Ic SNe)
    sne_Ibc = [sn for sn in sne_list if (sn.sn_type_class() == 'Ibc' \
                                    or  sn.sn_type_class() == 'Ib' \
                                    or  sn.sn_type_class() == 'Ic')]
    print('Found', len(sne_Ibc), 'SNe of type Ib, Ic or Ibc')

    # Ic list
    sne_Ic = [sn for sn in sne_list if sn.sn_type_class() == 'Ic']
    print('Found', len(sne_Ic), 'SNe of type Ic')

    # II list
    sne_II = [sn for sn in sne_list if sn.sn_type_class() == 'II']
    print('Found', len(sne_II), 'SNe of type II')

    # uses the given offset of the SN in the file to calculate the distance ratio
    offset_method = "read"

    # plots the Ia line (blue dot dash)
    plot_Ia, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ia]),\
             normalize(list(range(0, len(sne_Ia)))),\
             color = 'b', linestyle = '-.', linewidth = 2, label = 'Ia')

    # plots the Ib line (magenta dot dash)
    plot_Ib, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ib]),\
             normalize(list(range(0, len(sne_Ib)))),\
             color = 'm', linestyle = '-.', linewidth = 2, label = 'Ib')

    # plots the Ibc line (black dashed)
    plot_Ibc, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ibc]),\
             normalize(list(range(0, len(sne_Ibc)))),\
             color = 'k', linestyle = '--', linewidth = 2, label = 'Ibc')

    # plots the Ic line (gray solid)
    plot_Ic, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ic]),\
             normalize(list(range(0, len(sne_Ic)))),\
             color = '#777777', linestyle = '-', linewidth = 2, label = 'Ic')

    # plots the II line (red solid)
    plot_II, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_II]),\
             normalize(list(range(0, len(sne_II)))),\
             color = 'r', linestyle = '-', linewidth = 2, label = 'II')

    # x from 0.05 - 1.4 (distance ratio)
    plt.xlim(0.05, 1.4)

    # y from 0 - 1.1 ( cumulative frequency
    plt.ylim(0, 1.1)

    plt.legend(handles = [plot_II, plot_Ia, plot_Ib, plot_Ibc, plot_Ic], loc=4)

    plt.xlabel('Rsn / Rgal')
    plt.ylabel('Cumulative Fraction')
    plt.title('Duplicated Figure 8')

    # displays the plot
    plt.show()


def log_bins(low, high, n):
    return [low * 10 ** ((b / n) * math.log10(high/low)) for b in range(0, n+1)]


# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate(galaxy_bin, sn_type):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:
        if sn_type == 'Ia':
            total_sne = total_sne + g.sne_Ia                # stellar mass is in 10^10 Msun
            control_time_sum = control_time_sum + g.tc_Ia * g.stellar_mass_Lum * 10**10
        elif sn_type == 'SE':
            total_sne = total_sne + g.sne_SE
            control_time_sum = control_time_sum + g.tc_SE * g.stellar_mass_Lum * 10**10
        elif sn_type == 'II':
            total_sne = total_sne + g.sne_II
            control_time_sum = control_time_sum + g.tc_II * g.stellar_mass_Lum * 10**10
        else:
            return None

    return total_sne / control_time_sum


# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_outskirts(galaxy_bin, sn_type):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:
        num_Ia = 0
        num_SE = 0
        num_II = 0

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            curr_sn_type = sn.sn_type_class()

            # count the number of sne of a given type in the outskirts
            if sn.distance_ratio("read") >= 1.0:
                if curr_sn_type == 'Ia':
                    num_Ia = num_Ia + 1
                elif curr_sn_type == 'Ib' or curr_sn_type == 'Ic' or curr_sn_type == 'Ibc':
                    num_SE = num_SE + 1
                elif curr_sn_type == 'II':
                    num_II = num_II + 1

        if sn_type == 'Ia':
            total_sne = total_sne + num_Ia                  # stellar mass is in 10^10 Msun
            control_time_sum = control_time_sum + g.tc_Ia * g.stellar_mass_Lum * 10**10
        elif sn_type == 'SE':
            total_sne = total_sne + num_SE
            control_time_sum = control_time_sum + g.tc_SE * g.stellar_mass_Lum * 10**10
        elif sn_type == 'II':
            total_sne = total_sne + num_II
            control_time_sum = control_time_sum + g.tc_II * g.stellar_mass_Lum * 10**10
        else:
            return None


    return total_sne / control_time_sum


# returns the count of a certain type of supernova in a bin of galaxies
def sn_count(galaxy_bin, sn_type):
    total_sne = 0

    for g in galaxy_bin:
        if sn_type == 'Ia':
            total_sne = total_sne + g.sne_Ia                # stellar mass is in 10^10 Msun
        elif sn_type == 'SE':
            total_sne = total_sne + g.sne_SE
        elif sn_type == 'II':
            total_sne = total_sne + g.sne_II
        else:
            return None
    return total_sne


def bin_mean_mass(galaxy_bin):
    return sum([curr_gal.stellar_mass_Lum for curr_gal in galaxy_bin]) / len(galaxy_bin)


def calc_sne_rate():
    # parses the full sample of supernovae
    gal_dict = {}
    supernova_data.parse_galaxy_file(gal_dict, 'table2.dat')

    num_full_optimal = sum([1 for curr_gal in gal_dict.values() if curr_gal.full_optimal])
    num_not_optimal = sum([1 for curr_gal in gal_dict.values() if not curr_gal.full_optimal])

    print('full optimal', num_full_optimal)
    print('not optimal', num_not_optimal)

    # gets the number of each type of SNe (in the full-optimal sample only
    num_sne_Ia = sum([curr_gal.sne_Ia for curr_gal in gal_dict.values() if curr_gal.full_optimal])
    print('Samples SNe Ia:', num_sne_Ia)

    num_sne_SE = sum([curr_gal.sne_SE for curr_gal in gal_dict.values() if curr_gal.full_optimal])
    print('Samples SNe SE:', num_sne_SE)

    num_sne_II = sum([curr_gal.sne_II for curr_gal in gal_dict.values() if curr_gal.full_optimal])
    print('Samples SNe II:', num_sne_II)

    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal]

    n_bins = 1

    #                     10^8 Msun 10^12 Msun
    bin_limits = log_bins(0.01, 100, n_bins)

    galaxy_bins = []

    # for each mass bin
    for bin_num in range(0, n_bins):
        curr_bin_min = bin_limits[bin_num]
        curr_bin_max = bin_limits[bin_num + 1]

        galaxy_bin = []

        # for each galaxy
        for gal in list_galaxies:

            # check if the galaxy is in the right mass range
            if curr_bin_min < gal.stellar_mass_Lum and\
               gal.stellar_mass_Lum < curr_bin_max:
                galaxy_bin.append(gal)

        galaxy_bins.append(galaxy_bin)

    print(' *** SNe Ia *** ')

    sne_types = ['Ia', 'SE', 'II']

    for curr_sn_type in sne_types:

        mean_masses = [0] * n_bins
        sne_rates = [0] * n_bins

        print('Bin', 'low limit mass', 'high limit mass', 'sne')
        for bin_num in range(0, n_bins):
            #print('Bin: ' + str(bin_num),\
            #      round(bin_limits[bin_num], 3),\
            #      round(bin_limits[bin_num + 1], 3),\
            #      sne_bins[bin_num],\
            #      sn_rate(galaxy_bins[bin_num], curr_sn_type) * 10**12,\
            #      sep = '\t')

            mean_masses[bin_num] = bin_mean_mass(galaxy_bins[bin_num]) * 10**10
            sne_rates[bin_num] = sn_rate(galaxy_bins[bin_num], curr_sn_type) * 10**12

            print('count:', sn_count(galaxy_bins[bin_num], curr_sn_type))
            print('mass:', mean_masses[bin_num], end = '\t')
            print('rate:', sne_rates[bin_num])

        plt.scatter(mean_masses, sne_rates)
        plt.yscale('log')
        plt.xscale('log')

        plt.xlim(10**8, 10**12)

        plt.ylim(10**-2, 10**1)

        plt.xlabel('Solar Masses')
        plt.ylabel('SNuM')

        plt.title('SNe ' + curr_sn_type)

        plt.show()


def outskirts_rate_vs_mass():
    # parses the full sample of supernovae
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat')

    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    # break the galaxies into bins

    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal]

    n_bins = 1

    #                     10^8 Msun 10^12 Msun
    bin_limits = log_bins(0.01, 100, n_bins)

    galaxy_bins = []

    # for each mass bin
    for bin_num in range(0, n_bins):
        curr_bin_min = bin_limits[bin_num]
        curr_bin_max = bin_limits[bin_num + 1]

        galaxy_bin = []

        # for each galaxy
        for gal in list_galaxies:

            # check if the galaxy is in the right mass range
            if curr_bin_min < gal.stellar_mass_Lum and\
               gal.stellar_mass_Lum < curr_bin_max:
                galaxy_bin.append(gal)

        galaxy_bins.append(galaxy_bin)



    # find the rate for each bin

    sne_types = ['Ia', 'SE', 'II']

    for curr_sn_type in sne_types:
        print('***', curr_sn_type, '***')

        mean_masses = [0] * n_bins
        sne_rates = [0] * n_bins

        print('Bin', 'low limit mass', 'high limit mass', 'sne')
        for bin_num in range(0, n_bins):

            mean_masses[bin_num] = bin_mean_mass(galaxy_bins[bin_num]) * 10**10
            # finds the rate of SNe just in the outskirts
            sne_rates[bin_num] = sn_rate_outskirts(galaxy_bins[bin_num], curr_sn_type) * 10**12

            print('mass:', mean_masses[bin_num], end = '\t')
            print('rate:', sne_rates[bin_num])

        plt.scatter(mean_masses, sne_rates)
        plt.yscale('log')
        plt.xscale('log')

        #plt.xlim(10**8, 10**12)

        #plt.ylim(10**-2, 10**1)

        plt.title('SNe ' + curr_sn_type)

        plt.show()


def plot_num_sne_vs_radius():
    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full-optimal.txt')

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


def freedman_diaconis_nbins(arr):
    h = 2 * np.subtract(*np.percentile(arr, [75, 25])) / len(arr)**(1/3)
    return math.ceil((max(arr) - min(arr)) / h)


def count_types_in_outskirts():

    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full-optimal.txt')

    outskirt_sne = [sn for sn in sne_list if sn.distance_ratio('read') >= 1.0]

    num_Ia = len([sn for sn in outskirt_sne if sn.sn_type_class() == 'Ia'])

    num_SE = len([sn for sn in outskirt_sne if sn.sn_type_class() == 'Ib' or sn.sn_type_class() == 'Ic' or sn.sn_type_class() == 'Ibc'])

    num_II = len([sn for sn in outskirt_sne if sn.sn_type_class() == 'II'])

    print('SNe Ia in outskirts:', num_Ia)
    print('SE SNe in outskirts:', num_SE)
    print('SNe II in outskirts:', num_II)


def __main__():
    print(" *** sne_outskirts.py *** ")

    #plot_sne_cumulative()

    calc_sne_rate()

    #outskirts_rate_vs_mass()

    #plot_num_sne_vs_radius()


__main__()

