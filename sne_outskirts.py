# Brennan Dell
# Supernova Rate in Outskirts

import time
import math
import matplotlib.pyplot as plt
import numpy as np

import galaxy_data
import supernova_data
import color_profile
import spitzer_profile


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

def pair_galaxies_and_colors(gal_dict, color_dict):
    print('Pairing galaxies with color profiles...')
    profiled_galaxy_counter = 0

    for galaxy_name in gal_dict.keys():
        curr_FUV_name = to_FUV_name(galaxy_name)
        if curr_FUV_name in color_dict.keys():
            gal_dict[galaxy_name].color_profile = color_dict[curr_FUV_name]
            profiled_galaxy_counter = profiled_galaxy_counter + 1

    print(profiled_galaxy_counter, 'galaxies were given color profiles')

def pair_galaxies_and_spitzer_profiles(gal_dict, spitzer_dict):
    print('Pairing galaxies with spitzer profiles...')
    profiled_galaxy_counter = 0

    for galaxy_name in gal_dict.keys():
        curr_spitzer_name = to_FUV_name(galaxy_name)
        if curr_spitzer_name in spitzer_dict.keys():
            gal_dict[galaxy_name].band_36 = spitzer_dict[curr_spitzer_name][0]
            gal_dict[galaxy_name].band_45 = spitzer_dict[curr_spitzer_name][1]

            profiled_galaxy_counter = profiled_galaxy_counter + 1

    print(profiled_galaxy_counter, 'galaxies were given color profiles')


def test_axis_parsing():
    # parses galaxy data
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    test_galaxy = gal_dict[list(gal_dict.keys())[0]]
    print(test_galaxy.name, test_galaxy.maj_axis_min, test_galaxy.min_axis_min)

def test_galaxy_mass_function():
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the color profiles
    color_dict = {}
    color_profile.parse_color_profile_file(color_dict, 'nearby_galaxy_fuv_radius.txt')

    pair_galaxies_and_colors(gal_dict, color_dict)

    profiled_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.color_profile != None]

    test_galaxies = profiled_galaxies[0:10]
    for test_galaxy in test_galaxies:
        print(test_galaxy.name)

        test_galaxy.construct_mass_profile()
        test_mass_function = test_galaxy.mass_profile

        print('*** Mass Function ***')
        print('Radius', 'Mass', sep = '\t')
        for point in test_mass_function:
            print(point[0], point[1], sep = '\t')

        print('CMF(0.5) =', '%.3e' % test_galaxy.cumulative_mass_function(0.5))
        print('CMF(1.0) =', '%.3e' % test_galaxy.cumulative_mass_function(1.0))

        drs = [point[0] for point in test_mass_function]
        masses = [point[1] for point in test_mass_function]

        plt.plot(drs, masses)
        plt.show()

def find_mass_function_limits():
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the color profiles
    color_dict = {}
    color_profile.parse_color_profile_file(color_dict, 'nearby_galaxy_fuv_radius.txt')

    pair_galaxies_and_colors(gal_dict, color_dict)

    profiled_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.color_profile != None]

    print(time.time())

    for curr_gal in profiled_galaxies:
        curr_gal.construct_mass_profile()

    print(time.time())

    highest_radii = sorted([gal.mass_profile[-1][0] for gal in profiled_galaxies])

    mean_highest_radius = sum(highest_radii) / float(len(highest_radii))

    print("Mean highest radius:", mean_highest_radius)

    plt.hist(highest_radii, bins = 10, range = (1.0, 2.0), cumulative = -1, normed = True)
    plt.show()

    plt.hist(highest_radii, bins = 50, cumulative = -1, normed = True)
    plt.show()

def test_spitzer_parsing():
    # parses galaxy data
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    spitzer_dict = {}
    spitzer_profile.parse_spitzer_file(spitzer_dict, 'radial_profiles.csv')

    pair_galaxies_and_spitzer_profiles(gal_dict, spitzer_dict)

    galaxies_with_spitzer = [g for g in gal_dict.values() if g.band_36 is not None]

    print('*** Sample Galaxy ***')
    print('Name', galaxies_with_spitzer[0].name)
    print('Band 3.6')
    for k in galaxies_with_spitzer[0].band_36.keys():
        print('\t', k, galaxies_with_spitzer[0].band_36[k].aperture_flux)
    print('Band 4.5')
    for k in galaxies_with_spitzer[0].band_45.keys():
        print('\t', k, galaxies_with_spitzer[0].band_45[k].aperture_flux)

def count_colors():
    # parses the full sample of supernovae
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat')

    # parses the full sample of supernovae
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full.txt')

    # parses the color profiles
    color_dict = {}
    color_profile.parse_color_profile_file(color_dict, 'nearby_galaxy_fuv_radius.txt')

    pair_galaxies_and_colors(gal_dict, color_dict)

    pair_galaxies_and_sne(gal_dict, sne_list)

    profiled_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.color_profile != None]

    num_sne_with_profiles = sum([curr_gal.sne_Ia + curr_gal.sne_SE + curr_gal.sne_II\
                             for curr_gal in profiled_galaxies])

    sne_with_profiles = []
    for curr_gal in profiled_galaxies:
        sne_with_profiles = sne_with_profiles + curr_gal.supernovae

    print(len(sne_with_profiles), 'supernovae are in galaxies with color information')

    for curr_sn in sne_with_profiles:
        if curr_sn.distance_ratio("read") > 1.0:
            print(curr_sn.galaxy_name, curr_sn.sn_name, curr_sn.distance_ratio("read"), curr_sn.sn_type, sep = '\t')

def to_FUV_name(graur_name):
    if graur_name[0:3] == "ESO":
        if graur_name[0:6] == "ESO-LV":
            return None
        else:
            return "ESO" + graur_name[4:7] + "-" + graur_name[11:14]

    elif graur_name[0:2] == "IC":
        last_letter = ""
        if len(graur_name) > 7 and graur_name[7].isalpha():
            last_letter = graur_name[7]
        return "IC" + graur_name[3:7] + last_letter

    elif graur_name[0:3] == "NGC":
        last_letter = ""
        if ':' in graur_name:
            last_letter = ""    #ignore the one instance of NGC5457
        elif len(graur_name) > 8 and graur_name[8].isalpha():
            last_letter = graur_name[8]
        return "NGC" + graur_name[4:8] + last_letter

    elif graur_name[0:3] == "PGC":
        return "PGC" + graur_name[4:10]

    elif graur_name[0:3] == "UGC":
        if graur_name[0:4] == "UGCA":
            return None
        else:
            last_letter = ""
            if len(graur_name) > 9 and graur_name[9].isalpha():
                last_letter = graur_name[9]
            return "UGC" + graur_name[4:9] + last_letter
    else:
        return None


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

def sn_rate_general(sn_data):
    #sn_data is a list of tuples. Each tuple is: (Number SNe, mass of region (Msun), control time (yr))
    total_sne = 0
    control_time_sum = 0

    for g in sn_data:
        num_sne, stellar_mass, control_time = g

        total_sne = total_sne + num_sne
        control_time_sum = control_time_sum + control_time * stellar_mass

    if control_time_sum > 0:
        return total_sne / control_time_sum
    else:
        # return infinity if the mass of the calculated region is 0
        return float('inf')

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
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat')

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

def sn_rate_by_radius(num_bins, sn_type, low_limit, high_limit):
    # parses the full sample of supernovae
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    supernova_data.parse_sne_file(sne_list, 'sn-full.txt')

    # parses the color profiles
    color_dict = {}
    color_profile.parse_color_profile_file(color_dict, 'nearby_galaxy_fuv_radius.txt')

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

def __main__():
    print(" *** sne_outskirts.py *** ")

    #plot_sne_cumulative()

    #calc_sne_rate()

    #outskirts_rate_vs_mass()

    #plot_num_sne_vs_radius()

    #test_graur_to_FUV()

    #count_colors()

    #test_axis_parsing()

    #test_galaxy_mass_function()

    #find_mass_function_limits()

    #sn_rate_by_radius(2, 'Ia',0.0, 2.0)
    #sn_rate_by_radius(2, 'SE', 0.0, 2.0)
    #sn_rate_by_radius(2, 'II', 0.0, 2.0)

    test_spitzer_parsing()

__main__()

