import scipy.stats

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

#used for the rate by radius calculation.
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

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_total(galaxy_bin, sn_type):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:
        num_Ia = 0
        num_SE = 0
        num_II = 0

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            curr_sn_type = sn.sn_type_class()

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

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_outskirts_all_types(galaxy_bin, sn_type_throw_away):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            # count the number of sne of all types in the outskirts
            if sn.distance_ratio("read") >= 1.0:
                total_sne = total_sne + 1

        control_time_sum = control_time_sum + g.mean_control_time() * g.stellar_mass_Lum * 10**10

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_total_all_types(galaxy_bin):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            # count the number of sne of all types
            total_sne = total_sne + 1

        control_time_sum = control_time_sum + g.mean_control_time() * g.stellar_mass_Lum * 10 ** 10

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (
                sne_err_high - total_sne) / control_time_sum



def sn_rate_inside_all_types(galaxy_bin):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            # count the number of sne of all types in the outskirts
            if sn.distance_ratio("read") < 1.0:
                total_sne = total_sne + 1

        control_time_sum = control_time_sum + g.mean_control_time() * g.stellar_mass_Lum * 10 ** 10

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (
                sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_inside(galaxy_bin, sn_type):
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
            if sn.distance_ratio("read") < 1.0:
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

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum


# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_lum_outskirts_all_types(galaxy_bin):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            # count the number of sne of all types in the outskirts
            if sn.distance_ratio("read") >= 1.0:
                total_sne = total_sne + 1

        control_time_sum = control_time_sum + g.mean_control_time() * g.k_lum * 10**10

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_lum_total_all_types(galaxy_bin):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            # count the number of sne of all types
            total_sne = total_sne + 1

        control_time_sum = control_time_sum + g.mean_control_time() * g.k_lum * 10 ** 10

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (
                sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_lum_total(galaxy_bin, sn_type):
    total_sne = 0
    control_time_sum = 0
    for g in galaxy_bin:
        num_Ia = 0
        num_SE = 0
        num_II = 0

        # count the number of supernovae for this galaxy in the outskirts
        for sn in g.supernovae:
            curr_sn_type = sn.sn_type_class()

            if curr_sn_type == 'Ia':
                num_Ia = num_Ia + 1
            elif curr_sn_type == 'Ib' or curr_sn_type == 'Ic' or curr_sn_type == 'Ibc':
                num_SE = num_SE + 1
            elif curr_sn_type == 'II':
                num_II = num_II + 1

        if sn_type == 'Ia':
            total_sne = total_sne + num_Ia                  # stellar mass is in 10^10 Msun
            control_time_sum = control_time_sum + g.tc_Ia * g.k_lum * 10**10
        elif sn_type == 'SE':
            total_sne = total_sne + num_SE
            control_time_sum = control_time_sum + g.tc_SE * g.k_lum * 10**10
        elif sn_type == 'II':
            total_sne = total_sne + num_II
            control_time_sum = control_time_sum + g.tc_II * g.k_lum * 10**10
        else:
            return None

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum

# returns the rate of a certain type of supernova in a bin of galaxies
def sn_rate_lum_outskirts(galaxy_bin, sn_type):
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
            control_time_sum = control_time_sum + g.tc_Ia * g.k_lum * 10**10
        elif sn_type == 'SE':
            total_sne = total_sne + num_SE
            control_time_sum = control_time_sum + g.tc_SE * g.k_lum * 10**10
        elif sn_type == 'II':
            total_sne = total_sne + num_II
            control_time_sum = control_time_sum + g.tc_II * g.k_lum * 10**10
        else:
            return None

    sne_err_low, sne_err_high = scipy.stats.poisson.interval(0.68, total_sne)

    return total_sne / control_time_sum, (total_sne - sne_err_low) / control_time_sum, (sne_err_high - total_sne) / control_time_sum