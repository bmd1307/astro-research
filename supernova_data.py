import math

class Supernova:
    """Representation of a Supernova"""
    """Contains data from Table 4 from Leaman et al. (2010)"""
    """A Supernova object contains:"""
    """ sn_name                 - The supernova's name (ex. 1998S) """
    """ sn_offset_given         - The string representation of the offset of the SN from its host's nucleus in arcsecs (ex. '16.0W 4.5N') """
    """ sn_magnitude            - The supernova's magnitude """
    """ sn_ra_deg               - The supernova's right ascension in degrees (converted from hms in the file)"""
    """ sn_dec_deg              - The supernova's declination in degrees (converted from dms in the file)"""
    """ sn_type                 - The type of the supernova (Ia, Ib, ...)"""
    """ galaxy_name             - The name of the host galaxy (underscore format, ex. 'NGC_3877')"""
    """ galaxy_ra_deg           - The right ascension of the host galaxy in degrees (converted from hours format in the file)"""
    """ galaxy_dec_deg          - The declination of the host galaxy in degrees"""
    """ galaxy_dist_mpc         - The distance to the host galaxy in Mpc"""
    """ galaxy_hubble_type      - The hubble type of host galaxy (coded as an integer, see Leaman et al.)"""
    """ galaxy_maj_axis_min     - The major axis length of the host galaxy in arcmin"""
    """ galaxy_min_axis_min     - The minor axis length of the host galaxy in arcmin"""
    """ galaxy_incl_deg         - The inclination angle of the host galaxy"""
    """ galaxy_position_angle   - The posiiton angle of the host galaxy"""

    # Creates a Supernova object from a line from the supernova file
    def __init__(self, file_line):
        line_vals = file_line.split()

        self.sn_name = line_vals[0]

        self.offset_given = line_vals[6] + line_vals[7] + ' ' + line_vals[8] + line_vals[9]

        self.sn_magnitude = float(line_vals[10])

        sn_ra_hrs = float(line_vals[12])
        sn_ra_min = float(line_vals[13])
        sn_ra_sec = float(line_vals[14])
        # converts the hms format of the RA in the file to degrees
        self.sn_ra_deg = (sn_ra_hrs + (sn_ra_min / 60.0) + (sn_ra_sec / 3600.0)) * (360.0 / 24.0)

        # converts the dms format of the dec to degrees
        sn_dec_deg = float(line_vals[15])
        sn_dec_min = float(line_vals[16])
        sn_dec_sec = float(line_vals[17])
        self.sn_dec_deg = sn_dec_deg + math.copysign((sn_dec_min / 60.0) + (sn_dec_sec / 3600.0), sn_dec_deg)

        self.sn_type = line_vals[19]

        self.galaxy_name = line_vals[24]

        gal_ra = float(line_vals[25])
        # converts the hours format of the RA in the file to degrees
        self.galaxy_ra_deg = gal_ra * (360.0 / 24.0)

        self.galaxy_dec_deg = float(line_vals[26])

        self.galaxy_dist_mpc = float(line_vals[29])

        # the hubble type is coded as an integer (1-8)
        self.galaxy_hubble_type = int(line_vals[31])

        self.galaxy_maj_axis_min = float(line_vals[33])

        self.galaxy_min_axis_min = float(line_vals[34])

        self.galaxy_incl_deg = float(line_vals[35])

        self.galaxy_pa_deg = float(line_vals[37])

        self.galaxy_k_mag = float(line_vals[43]) if line_vals[43] != '99.999' else None
        self.galaxy_k_err = float(line_vals[44]) if self.galaxy_k_mag is not None else None

    # calculates the offset of the SN from the host galaxy by subtracting their coordinates
    # returns a string of the offset in arcsec (ex. '16.0W 4.5N')
    def offset(self):
        # get the ra and dec differences for the SN
        delta_ra_sec = (self.sn_ra_deg - self.galaxy_ra_deg) * 3600
        delta_dec_sec = (self.sn_dec_deg - self.galaxy_dec_deg) * 3600

        # round the ra to one decimal point of one arcsec
        offset_ra = str(abs(round(delta_ra_sec, 1)))
        # append W if the ra offset is negative, E if positive
        if delta_ra_sec < 0:
            offset_ra = offset_ra + 'W'
        else:
            offset_ra = offset_ra + 'E'

        # round the dec to one decimal point of one arcsec
        # append S if the ra offset is negative, N if positive
        offset_dec = str(abs(round(delta_dec_sec, 1)))
        if delta_dec_sec < 0:
            offset_dec = offset_dec + 'S'
        else:
            offset_dec = offset_dec + 'N'

        # assemble the parts and return the string
        return offset_ra + ' ' + offset_dec

    # calculate the distance ratio (Rsn / Rgal)
    # offset_type sets the calculation mode:
    #   "calc" - the offset of the SN is calculated by subtracting the SN ra and dec from the galaxy's ra and dec
    #   "read" - uses the offset as provided in the file (not always the same as the offset in "calc" mode)
    def distance_ratio(self, offset_type):
        axis_ratio = self.galaxy_maj_axis_min / self.galaxy_min_axis_min

        # gets the offset in arcsec
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

        # calculates the offset vector after it has been rotated so the galaxy's major axis is vertical, then dilates the x component by the axis ratio
        adjusted_offset_x = axis_ratio * \
                            (offset_x * cosd(self.galaxy_pa_deg) - \
                             offset_y * sind(self.galaxy_pa_deg))

        adjusted_offset_y = (offset_x * sind(self.galaxy_pa_deg) + \
                             offset_y * cosd(self.galaxy_pa_deg))

        maj_axis_sec = 60 * self.galaxy_maj_axis_min

        # returns the length of the adjusted offset vector (hypot(Ox, Oy)) divided by the galaxy's radius (1/2 the major axis)
        return math.hypot(adjusted_offset_x, adjusted_offset_y) / (0.5 * maj_axis_sec)

    # Determines the general SN class of this supernova
    # Ex. a SN type of 'IatXF' is converted to 'Ia'
    # The main classes of SNe are: '?', 'Ia', 'Ib', 'Ibc', 'Ic', and 'II'
    # anything after these main classes is ignored in the return string
    def sn_type_class(self):
        if self.sn_type[0] == '?':
            # if '?' in self.sn_type:
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


# a helper function: returns the sin of theta where theta is degrees (math.sin takes radians)
def sind(theta):
    return math.sin(math.radians(theta))


# a helper function: returns the cos of theta where theta is degrees (math.cos takes radians)
def cosd(theta):
    return math.cos(math.radians(theta))


# returns the pythagorean distance between two strings in offset format ('16.0W 4.5N')
def offset_distance(os1, os2):
    os1_c = offset_to_coords(os1)
    os2_c = offset_to_coords(os2)

    return math.hypot(os1_c[0] - os2_c[0], os1_c[1] - os2_c[1])


# returns a tuple of the difference between two offset strings ('16.0W 4.5N')
def offset_distance_manhattan(os1, os2):
    os1_c = offset_to_coords(os1)
    os2_c = offset_to_coords(os2)

    return (round(os1_c[0] - os2_c[0], 1), round(os1_c[1] - os2_c[1], 1))


# converts an offset string into a tuple in arcsec
# ex. '16.0W 4.5N' -> (-16.0, 4.5)
def offset_to_coords(os):
    os_vals = os.split()

    os_x = float(os_vals[0][:-1])
    if os_vals[0][-1] == 'W':
        os_x = -os_x

    os_y = float(os_vals[1][:-1])
    if os_vals[1][-1] == 'S':
        os_y = -os_y

    return (os_x, os_y)

# Reads each line in the file (file_name), creates a supernova object from it, and adds each supernova to the given array (sn_list)
def parse_sne_file(sn_list, file_name):
    print('Parsing file ' + file_name)

    sn_file = open(file_name)
    for line in sn_file:
        next_sn = Supernova(line)
        sn_list.append(next_sn)

    print('Parsed ' + str(len(sn_list)) + ' supernovae')