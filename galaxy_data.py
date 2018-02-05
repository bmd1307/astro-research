

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


def parse_galaxy_file(gal_dict, file_name):
    print('Parsing file ' + file_name)

    gal_file = open(file_name)
    for line in gal_file:
        next_gal = Galaxy(line)
        gal_dict[next_gal.name] = next_gal

    print('Parsed ' + str(len(gal_dict)) + ' galaxies')