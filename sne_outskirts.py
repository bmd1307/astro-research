# Brennan Dell
# Supernova Rate in Outskirts

import math
import matplotlib.pyplot as plt

class Supernova:
    """Representation of a Supernova"""
    def __init__(self, file_line):
        line_vals = file_line.split()
        
        self.sn_name            = line_vals[0]
        
        self.offset_given       = line_vals[6] + line_vals[7] + ' ' + line_vals[8] + line_vals[9]
        
        self.sn_magnitude       = float(line_vals[10]) 

        sn_ra_hrs               = float(line_vals[12])
        sn_ra_min               = float(line_vals[13])
        sn_ra_sec               = float(line_vals[14])
        # converts the hms format of the RA in the file to degrees
        self.sn_ra_deg          = (sn_ra_hrs + (sn_ra_min / 60.0) + (sn_ra_sec / 3600.0)) * (360.0 / 24.0)
        
        sn_dec_deg              = float(line_vals[15])
        sn_dec_min              = float(line_vals[16])
        sn_dec_sec              = float(line_vals[17])
        self.sn_dec_deg         = sn_dec_deg + math.copysign((sn_dec_min / 60.0) + (sn_dec_sec / 3600.0), sn_dec_deg)
        
        self.sn_type            = line_vals[19]

        self.galaxy_name        = line_vals[24]

        gal_ra                  = float(line_vals[25])
        # converts the hours format of the RA in the file to degrees
        self.galaxy_ra_deg      = gal_ra * (360.0 / 24.0)

        self.galaxy_dec_deg     = float(line_vals[26])

        self.galaxy_dist_mpc    = float(line_vals[29])
	
	# the hubble type is coded as an integer (1-8)
        self.galaxy_hubble_type = int(line_vals[31])

        self.galaxy_maj_axis_min = float(line_vals[33])

        self.galaxy_min_axis_min = float(line_vals[34])

        self.galaxy_incl_deg    = float(line_vals[35])

        self.galaxy_pa_deg      = float(line_vals[37])

        
    def offset(self):
        delta_ra_sec = (self.sn_ra_deg - self.galaxy_ra_deg) * 3600
        delta_dec_sec = (self.sn_dec_deg - self.galaxy_dec_deg) * 3600

        offset_ra = str(abs(round(delta_ra_sec, 1)))
        if delta_ra_sec < 0:
            offset_ra = offset_ra + 'W'
        else:
            offset_ra = offset_ra + 'E'

        offset_dec = str(abs(round(delta_dec_sec, 1)))
        if delta_dec_sec < 0:
            offset_dec = offset_dec + 'S'
        else:
            offset_dec = offset_dec + 'N'

        return offset_ra + ' ' + offset_dec

    def distance_ratio(self, offset_type):
        axis_ratio = self.galaxy_maj_axis_min / self.galaxy_min_axis_min

        offset_x = None
        offset_y = None
        if offset_type == "calc":
            offset_x = (self.sn_ra_deg - self.galaxy_ra_deg) * 3600
            offset_y = (self.sn_dec_deg - self.galaxy_dec_deg) * 3600
        elif offset_type == "read":
            offset_vals = offset_to_coords(self.offset_given)

            offset_x = offset_vals[0]
            offset_y = offset_vals[1]
        else:
            print('invalid offset type in distance_ratio')
        
        adjusted_offset_x = axis_ratio * \
                            (offset_x * cosd(self.galaxy_pa_deg) - \
                             offset_y * sind(self.galaxy_pa_deg))

        adjusted_offset_y = (offset_x * sind(self.galaxy_pa_deg) + \
                             offset_y * cosd(self.galaxy_pa_deg))

        maj_axis_sec = 60 * self.galaxy_maj_axis_min

        return math.hypot(adjusted_offset_x, adjusted_offset_y) / (0.5 * maj_axis_sec)

    def sn_type_class(self):
        if self.sn_type[0] == '?':
        #if '?' in self.sn_type:
            return '?'
        elif self.sn_type[:2] == 'II':
            return 'II'
        elif self.sn_type[1] == 'a':
            return 'Ia'
        elif self.sn_type[1] == 'c':
            return 'Ic'
        elif self.sn_type[1] == 'b':
            if len(self.sn_type) > 2 and self.sn_type[2] == 'c':
                return 'Ibc'
            else:
                return 'Ib'

def sind(theta):
    return math.sin(math.radians(theta))

def cosd(theta):
    return math.cos(math.radians(theta))

def offset_distance(os1, os2):
    os1_c = offset_to_coords(os1)
    os2_c = offset_to_coords(os2)

    return math.hypot(os1_c[0] - os2_c[0], os1_c[1] - os2_c[1])

def offset_distance_manhattan(os1, os2):
    os1_c = offset_to_coords(os1)
    os2_c = offset_to_coords(os2)

    return (round(os1_c[0] - os2_c[0], 1), round(os1_c[1] - os2_c[1], 1))

def offset_to_coords(os):
    os_vals = os.split()
        
    os_x = float(os_vals[0][:-1])
    if os_vals[0][-1] == 'W':
        os_x = -os_x

    os_y = float(os_vals[1][:-1])
    if os_vals[1][-1] == 'S':
        os_y = -os_y

    return (os_x, os_y)

def parse_sne_file(sn_list, file_name):
    print('Parsing file ' + file_name)

    sn_file = open(file_name)
    for line in sn_file:
        sn_list.append(Supernova(line))

    print('Parsed ' + str(len(sn_list)) + ' supernovae')

def normalize(arr):

    max_val = max(arr)
    min_val = min(arr)

    return [(n - min_val) / (max_val - min_val) for n in arr]

def __main__():
    print(" *** sne_outskirts.py *** ")

    sne_list = []

    parse_sne_file(sne_list, 'sn-full.txt')

    num_all_sne = len(sne_list)
    
    # remove the galaxies with a error in the position angle (-99.999)
    sne_list = [sn for sn in sne_list if sn.galaxy_pa_deg >= 0.0]

    print('Removed', str(num_all_sne - len(sne_list)), ' SNe with host galaxies that had bad position angles')

    sne_Ia = [sn for sn in sne_list if sn.sn_type_class() == 'Ia']
    print('Found', len(sne_Ia), 'SNe of type Ia')

    sne_Ib = [sn for sn in sne_list if sn.sn_type_class() == 'Ib']
    print('Found', len(sne_Ib), 'SNe of type Ib')

    # The Ibc plot in this graph is the Ibs, the Ics and Ibcs
    sne_Ibc = [sn for sn in sne_list if (sn.sn_type_class() == 'Ibc' \
                                    or  sn.sn_type_class() == 'Ib' \
                                    or  sn.sn_type_class() == 'Ic')]
    print('Found', len(sne_Ibc), 'SNe of type Ib, Ic or Ibc')

    sne_Ic = [sn for sn in sne_list if sn.sn_type_class() == 'Ic']
    print('Found', len(sne_Ic), 'SNe of type Ic')

    sne_II = [sn for sn in sne_list if sn.sn_type_class() == 'II']
    print('Found', len(sne_II), 'SNe of type II')

    offset_method = "read"

    plot_Ia, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ia]),\
             normalize(list(range(0, len(sne_Ia)))),\
             color = 'b', linestyle = '-.', linewidth = 2, label = 'Ia')

    plot_Ib, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ib]),\
             normalize(list(range(0, len(sne_Ib)))),\
             color = 'm', linestyle = '-.', linewidth = 2, label = 'Ib')

    plot_Ibc, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ibc]),\
             normalize(list(range(0, len(sne_Ibc)))),\
             color = 'k', linestyle = '--', linewidth = 2, label = 'Ibc')
    
    plot_Ic, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_Ic]),\
             normalize(list(range(0, len(sne_Ic)))),\
             color = '#777777', linestyle = '-', linewidth = 2, label = 'Ic')

    plot_II, = plt.step(sorted([sn.distance_ratio(offset_method) for sn in sne_II]),\
             normalize(list(range(0, len(sne_II)))),\
             color = 'r', linestyle = '-', linewidth = 2, label = 'II')

    plt.xlim(0.05, 1.4)
    plt.ylim(0, 1.1)

    plt.legend(handles = [plot_II, plot_Ia, plot_Ib, plot_Ibc, plot_Ic], loc=4)

    plt.xlabel('Rsn / Rgal')
    plt.ylabel('Cumulative Fraction')
    plt.title('Duplicated Figure 8')

    plt.show()
    
__main__()
    
