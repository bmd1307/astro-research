# Brennan Dell
# Supernova Rate in Outskirts

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
from experiments import *


# ******************************************************************************
# ***                          Helper Functions                              ***
# ******************************************************************************

# normalizes each value in a list to 0.0 - 1.0
# ex [1, 2, 3, 5] -> [0, 0.25, 0.5, 1.0]
# param: arr - the array to normalize
def normalize(arr):

    max_val = max(arr)
    min_val = min(arr)

    return [(n - min_val) / (max_val - min_val) for n in arr]


# produces an array of histogram bin limits increasing exponentially
# ex log_bins(1, 1000, 3) -> [1, 10, 100, 1000]
# param: low - the lower limit for the histogram
# param: high - the upper limit for the histogram
# param: n - the number of bins
def log_bins(low, high, n):
    return [low * 10 ** ((b / n) * math.log10(high/low)) for b in range(0, n+1)]


# returns the count of a certain type of supernova in a bin (list) of galaxies
# The valid types of supernovae: Ia, SE, and II
# param: galaxy_bin - a list of Galaxy objects
# param: sn_type, a string for the supernova type to be counted (valid values: "Ia", "SE", "II")
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


# returns the total number of supernova in a galaxy sample which occur beyond the optical radius
def count_total_sne_outskirts(galaxy_bin):
    to_return = 0

    for curr_gal in galaxy_bin:
        for curr_sn in curr_gal.supernovae:
            if curr_sn.distance_ratio("read") > 1.0:
                to_return = to_return + 1

    return to_return


# returns the mean of the stellar masses of a list of galaxy objects
# param: galaxy_bin - a list of Galaxy objects
def bin_mean_mass(galaxy_bin):
    return sum([curr_gal.stellar_mass_Lum for curr_gal in galaxy_bin]) / len(galaxy_bin)


# returns the ideal bin size for the values of 'arr' according to the freedman-diaconis equation:
#   bin_width = 2 * IQR(data) * len(data)^(-1/3)
# param: arr - the data to put into a bin
def freedman_diaconis_nbins(arr):
    h = 2 * np.subtract(*np.percentile(arr, [75, 25])) / len(arr)**(1/3)
    return math.ceil((max(arr) - min(arr)) / h)


# fits a power function to a series of x and y data using the linear regression on the logarithms of the x and y data
# results of the regression are printed, and the parameters of the fucntion are returned
# param: xs - the x values of the data (the independent variable)
# param: ys - the y values of the data (the dependent variable)
# returns: a tuple containing the coefficient and the exponent for the power function
def fit_exp(xs, ys):
    logxs = [math.log10(curr_x) for curr_x in xs]
    logys = [math.log10(curr_y) for curr_y in ys]

    results = scipy.stats.linregress(logxs, logys)

    print('Fitting a line to the log-log plot:')
    print('\tSlope:', results.slope, sep = '\t')
    print('\tIntercept:', results.intercept, sep = '\t')
    print('\tR Value:', results.rvalue, sep='\t')
    print('\tP Value:', results.pvalue, sep='\t')
    print('\tSTD Err:', results.stderr, sep='\t')

    print('\tCoefficient:', 10**results.intercept)
    print('\tExponent:', results.slope)

    return 10**results.intercept, results.slope


# ******************************************************************************
# ***                          Primary Functions                             ***
# ******************************************************************************


# Recreates the plot in Figure 8 of Leaman et al., the cumulative distributions of the SN types (Ia, Ib, Ibc, Ic, II)
def plot_sne_cumulative():

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

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


# Recreates the supernova rate calculations by leaman and graur
# Produces three plots of the supernova rates for type Ia, SE and II supernovae
# param: n_bins - the number of bins to group
def calc_sne_rate(n_bins):
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

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


# Plots a step function of the number of supernovae vs stellar mass
# A line is plotted for the total supernovae and the outskirts supernovae (beyond r25)
# param: n_bins - the number of bins (default 10)
def hist_total_sne_by_stellar_mass(n_bins = 10):
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    # limits the list of galaxies to the full_optimal
    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal]

    # constructs a list of the log of the stellar masses (if the mass is present) for galaxies with SNe
    log_masses = [math.log10(curr_gal.stellar_mass_Lum) for curr_gal in list_galaxies if curr_gal.stellar_mass_Lum > 0.0 and curr_gal.count_outskirts_sne() > 0]
    print('Optimal bin number:', freedman_diaconis_nbins(log_masses))

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
            if curr_bin_min < gal.stellar_mass_Lum and \
                    gal.stellar_mass_Lum < curr_bin_max:
                galaxy_bin.append(gal)

        galaxy_bins.append(galaxy_bin)

    sne_tot_bins = []
    sne_out_bins = []

    for gs in galaxy_bins:
        curr_total_sne = 0
        curr_total_out_sne = 0

        for curr_gal in gs:
            curr_total_sne = curr_total_sne + curr_gal.count_total_sne()

            curr_total_out_sne = curr_total_out_sne + curr_gal.count_outskirts_sne()
        sne_tot_bins.append(curr_total_sne)
        sne_out_bins.append(curr_total_out_sne)

    print('SNe bins:', sne_tot_bins)
    print('Total SNe found:', sum(sne_tot_bins))

    print('SNe (outskirts) bins:', sne_out_bins)
    print('Total SNe (outskirts) found:', sum(sne_out_bins))

    for i in range(0, len(sne_out_bins)):
        print('10^' + str(10.0 + math.log10(bin_limits[i])),\
              '10^' + str(10.0 + math.log10(bin_limits[i + 1])),\
              sne_tot_bins[i],\
              sne_out_bins[i], \
              '%.3g' % (1.0 * sne_out_bins[i] / sne_tot_bins[i]), \
              sep = '\t')

    plt.step(bin_limits[:-1], sne_tot_bins, where='post')
    plt.step(bin_limits[:-1], sne_out_bins, where='post')

    plt.legend(["Total SNe", 'Outskirts SNe'])

    #don't show the title for figures used in the paper, use captions instead
    #plt.title('Supernova Frequency vs Galaxy Stellar Mass')
    plt.xlabel('Stellar Mass ($10^{10} M_{\odot}$)')
    plt.ylabel('Number of SNe')
    plt.xscale('log')
    plt.show()


# Calculates the supernova rate in the outskirts. Recreates the plots as seen in Graur's paper.
# The rate is calculated with fixed bins and a sliding bins. The number of fixed bins is n_bins, and the sliding bin is
# calculated with the same width
# param: n_bins - the number of bins to calculate the rate
# param: save_graph - saves the graph directly to a png file. File name = out_rate_bins<n_bins>.png. False by default
# param: verbose - prints a lot of extra information about how the calculation is progressing and detailed results
#                   True by default
# param: show_graph - displays the graph and pauses the program. (default is True)
# param: rate_function - the function which performs the rate calculation valid functions: sn_rate_outskirts (default),
#                        sn_rate_total.
# param: title - the title for the graph. '(bins = <n_bins>' is added to this name. default: 'Supernova Rate in
#                   Outskirts'
# param: yrange - the lower and upper limits of the y axis (default: [0.001, 30])
def total_sn_rate_outskirts(n_bins, save_graph = False, verbose = True, show_graph = True,\
                            rate_function = sn_rate_outskirts, title = 'Supernova Rate in Outskirts', yrange = [0.001, 30]):
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    # limits the list of galaxies to the full_optimal
    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal]

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
            if curr_bin_min < gal.stellar_mass_Lum and \
                    gal.stellar_mass_Lum < curr_bin_max:
                galaxy_bin.append(gal)

        galaxy_bins.append(galaxy_bin)

    rate_bins = []
    low_err_bins = []
    high_err_bins = []

    for curr_bin in galaxy_bins:
        curr_Ia_rate = rate_function(curr_bin, 'Ia')
        curr_SE_rate = rate_function(curr_bin, 'SE')
        curr_II_rate = rate_function(curr_bin, 'II')

        total_rate = curr_Ia_rate[0] + \
                     curr_SE_rate[0] + \
                     curr_II_rate[0]

        total_low = curr_Ia_rate[1] + \
                    curr_SE_rate[1] + \
                    curr_II_rate[1]

        total_high = curr_Ia_rate[2] + \
                     curr_SE_rate[2] + \
                     curr_II_rate[2]

        # convert to SNuM
        rate_bins.append(total_rate * 10**12)
        low_err_bins.append(total_low * 10**12)
        high_err_bins.append(total_high * 10**12)

    mean_masses = [bin_mean_mass(curr_bin) for curr_bin in galaxy_bins]

    # delete all the points with a rate of zero
    bad_indices = [n for n in range(0, len(rate_bins)) if rate_bins[n] == 0.0]

    mean_masses = [mean_masses[i] for i in range(0, len(mean_masses)) if i not in bad_indices]
    rate_bins = [rate_bins[i] for i in range(0, len(rate_bins)) if i not in bad_indices]
    low_err_bins = [low_err_bins[i] for i in range(0, len(low_err_bins)) if i not in bad_indices]
    high_err_bins = [high_err_bins[i] for i in range(0, len(high_err_bins)) if i not in bad_indices]

    if verbose:
        for i in range(0, len(rate_bins)):
            print('10^' + str(10.0 + math.log10(mean_masses[i])),\
                rate_bins[i],\
                sep = '\t')

    #plt.scatter(mean_masses, rate_bins, c = 'r')
    plt.errorbar(mean_masses, rate_bins, yerr = [high_err_bins, low_err_bins], color = 'r', zorder = 10)

    fit_coefficient, fit_exponent = fit_exp(mean_masses, rate_bins)

    rate = lambda stellar_mass_10: fit_coefficient * stellar_mass_10 ** fit_exponent

    fit_xs = [10**n for n in range(-2, 3)]
    fit_ys = [rate(x) for x in fit_xs]
    plt.plot(fit_xs, fit_ys, color = 'lime', zorder = 20)

    fit_func_str = ('%1.3f' % fit_coefficient) + ' * M* ^ ' + ('%1.3f' % fit_exponent)

    print('Fit Function (10^10 Msun -> SNuM):', fit_func_str)

    if verbose:
        print('Sliding bin')
    # add a plot for a sliding bin
    sliding_bin_low_limits = log_bins(10**-3, 10**2, 350)
    sliding_bin_width = (10**12 / 10**8) ** (1 / n_bins)
    mean_masses = []
    rate_bins = []
    low_err_bins = []
    high_err_bins = []

    for low_limit in sliding_bin_low_limits:
        high_limit = sliding_bin_width * low_limit

        curr_bin = [curr_gal for curr_gal in list_galaxies if low_limit < curr_gal.stellar_mass_Lum < high_limit]

        if len(curr_bin) == 0:
            continue

        curr_Ia_rate = rate_function(curr_bin, 'Ia')
        curr_SE_rate = rate_function(curr_bin, 'SE')
        curr_II_rate = rate_function(curr_bin, 'II')

        total_rate = curr_Ia_rate[0] + \
                     curr_SE_rate[0] + \
                     curr_II_rate[0]

        total_low = curr_Ia_rate[1] + \
                    curr_SE_rate[1] + \
                    curr_II_rate[1]

        total_high = curr_Ia_rate[2] + \
                     curr_SE_rate[2] + \
                     curr_II_rate[2]

        mean_masses.append(bin_mean_mass(curr_bin))
        rate_bins.append(total_rate * 10**12)
        low_err_bins.append(total_low * 10**12)
        high_err_bins.append(total_high * 10**12)


    # delete all the points with a rate of zero
    bad_indices = [n for n in range(0, len(rate_bins)) if rate_bins[n] == 0.0]

    mean_masses = [mean_masses[i] for i in range(0, len(mean_masses)) if i not in bad_indices]
    rate_bins = [rate_bins[i] for i in range(0, len(rate_bins)) if i not in bad_indices]
    low_err_bins = [low_err_bins[i] for i in range(0, len(low_err_bins)) if i not in bad_indices]
    high_err_bins = [high_err_bins[i] for i in range(0, len(high_err_bins)) if i not in bad_indices]

    if verbose:
        print('Number of bins for moving mass:', len(rate_bins))
    if verbose:
        for i in range(0, len(rate_bins)):
            print('Mean mass:', '%.3e' % mean_masses[i], '\tRate:', rate_bins[i])

    plt.plot(mean_masses, rate_bins, zorder = 5)

    upper_line = [rate_bins[i] + high_err_bins[i] for i in range(0, len(rate_bins))]
    lower_line = [rate_bins[i] - low_err_bins[i] for i in range(0, len(rate_bins))]

    plt.fill_between(mean_masses, lower_line, upper_line, facecolor = 'gainsboro', zorder = 1)

    #don't show this title, use caption instead in the paper
    #plt.title(title + ' (bins = ' + str(n_bins) + ')')
    leg = plt.legend(['Linear Fit to Fixed Bins\n(' + fit_func_str + ')', 'Sliding Bin',\
                'Sliding bin error (68% Poisson)',\
                'Fixed Bin'])
    leg.set_zorder(30)
    plt.xlabel('Stellar Mass ($10^{10} M_{\odot}$)')
    plt.ylabel('SNum ($10^{-12} M_{\odot}^{-1} yr^{-1}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([10**-2.5, 10**2.5])
    plt.ylim(yrange)

    if show_graph:
        plt.show()

    if save_graph:
        plt.savefig('..\\..\\rate_images\\out_rate_bins' + str(n_bins) + '.png', bbox_inches = 'tight')
        plt.gcf().clear()


# calculates the supernova rate in the dwarf galaxies in the galaxy sample
# prints the details of the calculation
def total_rate_dwarfs():
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    # limits the list of galaxies to the full_optimal
    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal and 0.0 < curr_gal.stellar_mass_Lum < 0.1]

    num_sne = sum([len(curr_gal.supernovae) for curr_gal in list_galaxies])

    mean_mass = bin_mean_mass(list_galaxies)

    lowest_mass = min([curr_gal.stellar_mass_Lum for curr_gal in list_galaxies])
    highest_mass = max([curr_gal.stellar_mass_Lum for curr_gal in list_galaxies])

    print('Found', len(list_galaxies), 'dwarf galaxies (M* < 10^9 Msun)')
    print('These galaxies host', num_sne, 'supernovae')
    print('Mean galaxy mass:', '%.3e' % (mean_mass * 1e10), 'Msun')
    print('Minimum galaxy mass:', '%.3e' % (lowest_mass * 1e10), 'Msun')
    print('Maximum galaxy mass:', '%.3e' % (highest_mass * 1e10), 'Msun')


    print('Hubble types:')

    print('\tElliptical:', len([0 for curr_gal in list_galaxies if curr_gal.hubble_type[0] == 'E']), sep = '\t')
    print('\tSpiral:', len([0 for curr_gal in list_galaxies if curr_gal.hubble_type[0] == 'S']), sep = '\t')
    print('\tIrregular:', len([0 for curr_gal in list_galaxies if curr_gal.hubble_type[0] == 'I']), sep = '\t')


    rate_Ia, low_Ia, high_Ia = sn_rate_total(list_galaxies, 'Ia')
    rate_SE, low_SE, high_SE = sn_rate_total(list_galaxies, 'SE')
    rate_II, low_II, high_II = sn_rate_total(list_galaxies, 'II')

    print('Type Ia rate:', '%1.3f' % (rate_Ia * 1e12), '(-' '%1.3f' % (low_Ia * 1e12), ', +' + \
          '%1.3f' % (high_Ia * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print('Type SE rate:', '%1.3f' % (rate_SE * 1e12), '(-' '%1.3f' % (low_SE * 1e12), ', +' + \
          '%1.3f' % (high_SE * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print('Type II rate:', '%1.3f' % (rate_II * 1e12), '(-' '%1.3f' % (low_II * 1e12), ', +' + \
          '%1.3f' % (high_II * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print()

    total_rate = rate_Ia + rate_SE + rate_II
    total_low = low_Ia + low_SE + low_II
    total_high = high_Ia + high_SE + high_II

    print('Total rate:', '%1.3f' % (total_rate * 1e12), '(-' '%1.3f' % (total_low * 1e12), ', +' + \
          '%1.3f' % (total_high * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')

    # calculates the number of SNe per millenium
    sne_per_mill = total_rate * mean_mass * 1e10 * 1e3
    sne_per_mill_low = total_low * mean_mass * 1e10 * 1e3
    sne_per_mill_high = total_high * mean_mass * 1e10 * 1e3

    print('SNe / millenium (rate * avg mass):', \
          '%1.3f' % sne_per_mill, '(-', '%1.3f' % sne_per_mill_low, ',+', '%1.3f' % sne_per_mill_high, ')',
          'SNe / (1000 yr) ')


# calculates the supernova rate in the outskirts of the spiral galaxies in the galaxy sample
# prints the details of the calculation
def total_rate_outskirts_spirals():
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae
    sne_list = []
    parse_sne_file(sne_list, 'sn-full.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    # limits the list of galaxies to the full_optimal with spirals larger than 10^9 Msun
    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.full_optimal and curr_gal.hubble_type[0] == 'S' and curr_gal.stellar_mass_Lum > 0.1]

    # gets mass data for the galaxies
    log_gal_masses = [math.log10(10**10 * curr_gal.stellar_mass_Lum) for curr_gal in list_galaxies]

    plt.hist(log_gal_masses, bins = 20)
    plt.title("Number of Spiral Galaxies vs Stellar Mass")
    plt.xlabel("Stellar Mass (log(Msun))")
    plt.ylabel("Number of Galaxies")
    plt.show()

    num_sne = count_total_sne_outskirts(list_galaxies)

    mean_mass = bin_mean_mass(list_galaxies)
    lowest_mass = min([curr_gal.stellar_mass_Lum for curr_gal in list_galaxies])
    highest_mass = max([curr_gal.stellar_mass_Lum for curr_gal in list_galaxies])

    print('Found', len(list_galaxies), 'spiral galaxies')
    print('These galaxies host', num_sne, 'supernovae')
    print('Mean galaxy mass:', '%.3e' % (mean_mass * 1e10), 'Msun')
    print('Minimum galaxy mass:', '%.3e' % (lowest_mass * 1e10), 'Msun')
    print('Maximum galaxy mass:', '%.3e' % (highest_mass * 1e10), 'Msun')

    rate_Ia, low_Ia, high_Ia = sn_rate_outskirts(list_galaxies, 'Ia')
    rate_SE, low_SE, high_SE = sn_rate_outskirts(list_galaxies, 'SE')
    rate_II, low_II, high_II = sn_rate_outskirts(list_galaxies, 'II')

    print('Type Ia rate:', '%1.3f' % (rate_Ia * 1e12), '(-' '%1.3f' % (low_Ia * 1e12), ', +' + \
          '%1.3f' % (high_Ia * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print('Type SE rate:', '%1.3f' % (rate_SE * 1e12), '(-' '%1.3f' % (low_SE * 1e12), ', +' + \
          '%1.3f' % (high_SE * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print('Type II rate:', '%1.3f' % (rate_II * 1e12), '(-' '%1.3f' % (low_II * 1e12), ', +' + \
          '%1.3f' % (high_II * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')
    print()

    total_rate = rate_Ia + rate_SE + rate_II
    total_low = low_Ia + low_SE + low_II
    total_high = high_Ia + high_SE + high_II

    print('Total rate:', '%1.3f' % (total_rate * 1e12), '(-' '%1.3f' % (total_low * 1e12), ', +' + \
          '%1.3f' % (total_high * 1e12), ')', '* 10^-12 SNe yr^-1 Msun^-1')

    # calculates the number of SNe per millenium
    sne_per_mill = total_rate * mean_mass * 1e10 * 1e3
    sne_per_mill_low = total_low * mean_mass * 1e10 * 1e3
    sne_per_mill_high = total_high * mean_mass * 1e10 * 1e3

    print('SNe / millenium (rate * avg mass):', \
          '%1.3f' % sne_per_mill, '(-', '%1.3f' % sne_per_mill_low, ',+', '%1.3f' % sne_per_mill_high, ')',
          'SNe / (1000 yr) ')


# runs total_rate_dwarfs and total_rate_outskirts_spirals
def compare_outskirts_to_dwarfs():
    print('*' * 80)
    print('Calculating the total supernova rate in dwarf galaxies')
    print('*' * 80)
    print()

    total_rate_dwarfs()

    print()
    print('*' * 80)
    print('Calculating the outskirts supernova rate in spiral galaxies')
    print('*' * 80)

    total_rate_outskirts_spirals()


# finds the radial dependence of supernovae. Groups the supernovae into an ideal number of bins by radius, and fits a
# power function to the supernova density as a function of radius. Displays a graph of the supernova rates
def sne_radial_data():
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    parse_sne_file(sne_list, 'sn-full-optimal.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    distance_ratios = [sn.distance_ratio('read') for sn in sne_list]

    distance_ratios = list(filter(lambda x: 0.5 < x < 4.0, distance_ratios))

    print(len(distance_ratios))

    print(freedman_diaconis_nbins(distance_ratios))

    counts, bin_lims =  np.histogram(distance_ratios, freedman_diaconis_nbins(distance_ratios))
    print('counts:', counts)
    print('bin_lims:', bin_lims)

    bin_centers = [0.5 * (bin_lims[i] + bin_lims[i + 1]) for i in range(0, len(bin_lims) - 1)]

    bin_width = (bin_lims[-1] - bin_lims[0]) / (len(bin_lims) - 1)

    good_indices = [i for i in range(0, len(counts)) if counts[i] != 0]

    counts = [counts[i] for i in good_indices]

    bin_centers = [bin_centers[i] for i in good_indices]

    # number of supernovae / (4 pi r dr). 4 pi r dr is the area of a ring with central radius 'r' and width '2dr'
    sn_densities = [counts[i] / (4 * math.pi * bin_centers[i] * (0.5 * bin_width)) for i in range(0, len(bin_centers))]

    print(sn_densities)

    print('*** Fitting a SN count to an exp curve ***')
    fit_exp(bin_centers, counts)

    print('*** Fitting a SN densities to an exp curve ***')
    density_coefficient, density_exponent = fit_exp(bin_centers, sn_densities)

    density_func = lambda r: density_coefficient * r ** density_exponent

    fit_func_str = ('%1.3f' % density_coefficient) + ' * (r / R25) ^ ' + ('%1.3f' % density_exponent)

    # xs are the distance ratios from 0 - 5
    fit_xs = [0.1 * n for n in range(1, 50)]
    fit_ys = [density_func(x) for x in fit_xs]
    plt.plot(fit_xs, fit_ys, color='lime', zorder=20)

    plt.title('Supernova Density vs Distance Ratio')

    plt.xlabel('Distance Ratio')
    plt.ylabel('Supernova Density')
    plt.xlim([0.4, 5.0])
    plt.ylim([0.1, 1000])
    plt.xscale('log')
    plt.yscale('log')
    #plt.scatter(bin_centers, counts, c='blue')
    plt.scatter(bin_centers, sn_densities, c='red')

    plt.legend(['Best fit: ' + fit_func_str, 'Supernova Density Bins'])

    plt.show()


# plots a histogram of the supernova frequency vs radius.
# the number of bins is chosen by the freedman-diaconis equation
def sne_radial_histogram():
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # parses the full sample of supernovae (contains those of '?' type)
    sne_list = []
    parse_sne_file(sne_list, 'sn-full-optimal.txt')

    pair_galaxies_and_sne(gal_dict, sne_list)

    distance_ratios = [sn.distance_ratio('read') for sn in sne_list]

    num_in_outskirts = count_total_sne_outskirts(gal_dict.values())

    print("Total SNe beyond R25:", num_in_outskirts)
    print("Total number of supernovae:", len(distance_ratios))

    plt.hist(distance_ratios, bins = freedman_diaconis_nbins(distance_ratios))
    # don't use title in the figures for the paper, caption instead
    #plt.title('Number of Supernovae vs Radius (r / $R_{25}$)')
    plt.xlabel('Galactocentric Radius (r / $R_{25}$)')
    plt.ylabel('Number of SNe')
    plt.show()


def galaxy_mass_histogram(log_scale = False):
    # parses the full sample of supernovae
    gal_dict = {}
    parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    # gets all the galaxies with a non zero mass
    list_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.stellar_mass_Lum > 0.0]

    # creates a list of stellar masses, converts
    galaxy_masses = [curr_gal.stellar_mass_Lum * 10**10 for curr_gal in list_galaxies]

    if log_scale:
        galaxy_masses = [math.log10(curr_mass) for curr_mass in galaxy_masses]

    n_bins = freedman_diaconis_nbins(galaxy_masses)

    print('Freedman Diaconis bins:', n_bins)

    plt.hist(galaxy_masses, bins=n_bins)

    plt.title('Stellar Mass Function of the Galaxy Sample')

    plt.xlabel('Log10(M*)')

    plt.ylabel('Number of Galaxies')

    plt.yscale('log')

    plt.show()

# call the desired function from this __main__ function
def __main__():
    print(" *** sne_outskirts.py *** ")

    #hist_total_sne_by_stellar_mass(10)

    #for n in range(2, 25):
    #    print('Saving image for n =', n)
    #    total_sn_rate_outskirts(n, save_graph=True, verbose=False, show_graph=False)

    #total_sn_rate_outskirts(10, rate_function=sn_rate_total, title = 'Total Supernova Rate vs Stellar Mass', yrange=[0.01, 300])
    #total_sn_rate_outskirts(10) # calculates outskirts by default

    #test_outskirts_mass_diaz_garcia()

    #sne_radial_data()

    #sne_radial_histogram()

    compare_outskirts_to_dwarfs()

    #galaxy_mass_histogram(log_scale=True)


__main__()

