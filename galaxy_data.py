

class Galaxy:
    # Creates a Supernova object from a line from the supernova file
    def __init__(self, file_line):
        self.name                   = file_line[0:21].strip()
        self.full_optimal           = bool(int(file_line[22]))
        self.sne_Ia                 = int(file_line[24])
        self.sne_SE                 = int(file_line[26])
        self.sne_II                 = int(file_line[28])
        self.tc_Ia                  = float(file_line[30:36])
        self.tc_SE                  = float(file_line[37:42])
        self.tc_II                  = float(file_line[43:49])
        self.stellar_mass_Lum       = float(file_line[50:58])
        self.distance               = float(file_line[76:82])
        self.hubble_type            = file_line[83:86]
        self.ra_deg                 = float(file_line[87:95])
        self.dec_deg                = float(file_line[96:104])

        # # These attributes don't exist for each galaxy in the .dat file
        # self.stellar_mass_galspec   = float(file_line[128:135])
        # self.sfr_SDSS               = float(file_line[136:143])
        # self.sfr_galspec            = float(file_line[144:153])
        # self.metallicity            = float(file_line[154:160])
        # self.sfr_petrosian          = float(file_line[161:168])
        # self.sfr_sersic             = float(file_line[169:176])

        self.supernovae             = None # initialize this to none to indicate that it hasn't been parsed yet, not that it has zero galaxies

        self.color_profile          = None
        self.mass_profile           = None

        self.maj_axis_min           = None # not present in the graur data
        self.min_axis_min           = None # these come from the Graur data

    def construct_mass_profile(self):
        if self.color_profile == None:
            raise AttributeError('A mass profile cannot be constructed for galaxies with no color profiles')

def parse_galaxy_file(gal_dict, main_file_name, axes_file_name):
    print('Parsing file ' + main_file_name)

    gal_file = open(main_file_name)
    for line in gal_file:
        next_gal = Galaxy(line)
        gal_dict[next_gal.name] = next_gal

    print('Parsed ' + str(len(gal_dict)) + ' galaxies')

    count_galaxies_no_axes = 0
    count_axes_data = 0

    axes_file = open(axes_file_name)
    for line in axes_file:
        curr_gal_name = line[10:33].strip()

        if curr_gal_name in gal_dict.keys():
            count_axes_data = count_axes_data + 1
            gal_dict[curr_gal_name].maj_axis_min = float(line[91:96])
            gal_dict[curr_gal_name].min_axis_min = float(line[99:104])
        else:
            count_galaxies_no_axes = count_galaxies_no_axes + 1

    print('Parsed axis data for', count_axes_data, 'galaxies.')

    if count_galaxies_no_axes > 0:
        print(count_galaxies_no_axes, 'were not given axis data')
    else:
        print('All galaxies have been given axis data')