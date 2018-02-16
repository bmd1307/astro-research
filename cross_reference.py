

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