import time
import math
import matplotlib.pyplot as plt
import numpy as np

import galaxy_data
import supernova_data
import color_profile
import spitzer_profile

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

def test_spitzer_profiles():
    gal_dict = {}
    galaxy_data.parse_galaxy_file(gal_dict, 'table2.dat', 'galaxy-full.txt')

    spitzer_dict = {}
    spitzer_profile.parse_spitzer_file(spitzer_dict, 'radial_profiles.csv')

    pair_galaxies_and_spitzer_profiles(gal_dict, spitzer_dict)

    profiled_galaxies = [curr_gal for curr_gal in gal_dict.values() if curr_gal.band_36 is not None]

    test_galaxies = profiled_galaxies[0:10]
    for test_galaxy in test_galaxies:
        print(test_galaxy.name)

        test_galaxy.construct_spitzer_profile()
        test_mass_function = test_galaxy.mass_profile

        print('*** Mass Function ***')
        print('Radius', 'Mass', sep='\t')
        for point in test_mass_function:
            print('%.3f' % point[0], '%.3e' % point[1], sep='\t')

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