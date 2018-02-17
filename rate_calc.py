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


    return total_sne / control_time_sum