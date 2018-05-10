import math

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

        # These values may not exist, so set them as none
        self.stellar_mass_galspec = None
        self.sfr_SDSS             = None
        self.sfr_galspec          = None
        self.metallicity          = None
        self.sfr_petrosian        = None
        self.sfr_sersic           = None

        self.k_mag                = None
        self.k_err                = None

        # # These attributes don't exist for each galaxy in the .dat file
        # so check that the value exists before writing it
        if(len(file_line) >= 135 and file_line[128:135].strip() != ''):
            self.stellar_mass_galspec   = float(file_line[128:135])
        if(len(file_line) >= 143 and file_line[136:143].strip() != ''):
            self.sfr_SDSS               = float(file_line[136:143])
        if(len(file_line) >= 153 and file_line[144:153].strip() != ''):
            self.sfr_galspec            = float(file_line[144:153])
        if(len(file_line) >= 160 and file_line[154:160].strip() != ''):
            self.metallicity            = float(file_line[154:160])
        if(len(file_line) >= 168 and file_line[161:168].strip() != ''):
            self.sfr_petrosian          = float(file_line[161:168])
        if(len(file_line) >= 176 and file_line[169:176].strip() != ''):
            self.sfr_sersic             = float(file_line[169:176])

        self.supernovae             = None # initialize this to none to indicate that it hasn't been parsed yet, not that it has zero galaxies

        self.color_profile          = None
        self.mass_profile           = None

        self.flux_36_profile        = None
        self.flux_45_profile        = None

        self.maj_axis_min           = None # not present in the graur data
        self.min_axis_min           = None # these come from the Graur data

        self.band_36                = None # band 36 data from the spitzer images
        self.band_45                = None  # band 36 data from the spitzer images

    def r25_pc(self):
        if self.maj_axis_min is None:
            return None
        return (0.5 * 60 * self.maj_axis_min) * (self.distance * 1e6) / 206265.0

    def construct_mass_profile(self):
        if self.color_profile == None:
            raise AttributeError('A mass profile cannot be constructed for galaxies with no color profiles')

        if self.maj_axis_min == None:
            raise AttributeError('A mass profile cannot be constructed for galaxies with no major axis measurement')

        self.mass_profile = []

        # get the indices of the color profile which contain IR data and sort them in ascending order
        indices = sorted([index for index in self.color_profile.keys() if self.color_profile[index].mu_36 != None])

        arcsec_to_distance_ratio = lambda sec: float(sec) / (self.maj_axis_min * 60) #convert the minutes to seconds

                                                    # convert Mpc to pc
        arcsec_to_parsec = lambda sec: float(sec) * (self.distance * 1000000.0) / 206265.0

        # accumulates the mass
        cumulative_mass = 0.0

        # assume all the mass is contained from the outermost ring inwards
        self.mass_profile.append((arcsec_to_distance_ratio(3), cumulative_mass))

        for rad_index in indices:
            curr_ring_radius = arcsec_to_parsec(rad_index)
            curr_ring_half_width = arcsec_to_parsec(3)

            curr_ellipticity = self.color_profile[rad_index].ellipticity

            # the area of a ring-shaped ellipse = (4pi*central radius * half_ring_width) * (1 - ellipticity)
            curr_ring_area_pc = 4 * math.pi * curr_ring_radius * curr_ring_half_width * (1 - curr_ellipticity)

            # calculates the mass in the ring in solar masses (surface density times area in square parsecs)
            curr_ring_mass = self.color_profile[rad_index].ring_mass_density() * curr_ring_area_pc

            # add the next point on the surface profile
            self.mass_profile.append((arcsec_to_distance_ratio(rad_index + 3), cumulative_mass + curr_ring_mass))

            # the new cumulative mass is the current cumulative mass - the current ring's mass
            cumulative_mass = cumulative_mass + curr_ring_mass

    def cumulative_mass_function(self, distance_ratio):
        if self.mass_profile is None:
            raise AttributeError("Cannot access a galaxy's CMF which hasn't been constructed.")

        drs = [point[0] for point in self.mass_profile]

        right_point_index = 0

        # find the index of the point in the mass function with the first distance ratio which is greater than the
        # requested distance ratio
        while right_point_index < len(drs) and drs[right_point_index] < distance_ratio:
            right_point_index = right_point_index + 1

        # the distance ratio is within the innermost ring (outside the domain of the CMF)
        if right_point_index == 0:
            return 0
        # if the distance ratio is beyond the outermost ring
        elif right_point_index == len(drs):
            return self.mass_profile[-1][1] # return the mass of the last point on the CMF
        # if the distance ratio is between two points on the CMF
        else:
            # assume the point at distance_ratio is on the line between its two nearest points
            left_dr = self.mass_profile[right_point_index - 1][0]
            right_dr = self.mass_profile[right_point_index][0]

            left_mass = self.mass_profile[right_point_index - 1][1]
            right_mass = self.mass_profile[right_point_index][1]

            slope = (right_mass - left_mass) / (right_dr - left_dr)

            return (distance_ratio - left_dr) * slope + left_mass

    def construct_spitzer_profile(self):
        if self.band_36 is None or self.band_45 is None:
            raise AttributeError('A spitzer mass profile cannot be constructed for galaxies with no IR band data')

        if self.maj_axis_min is None:
            raise AttributeError('A mass profile cannot be constructed for galaxies with no major axis measurement')

        arcsec_to_distance_ratio = lambda sec: float(sec) / (self.maj_axis_min * 60)  # convert the minutes to seconds

        # fluxes are in Jy, distance is in Mpc
        flux_to_mstar = lambda f36_Jy, f45_Jy, d_mpc: (10**5.65) * (f36_Jy**2.85) * (f45_Jy**-1.85) * (d_mpc / 0.05)**2

        common_indices = [index for index in self.band_36 if index in self.band_45]

        indices_with_data = sorted([index for index in common_indices\
                                if self.band_36[index] is not None and self.band_45[index] is not None])

        # the mass at distance ratio 0 is 0
        self.mass_profile = [(0.0, 0.0)]

        for index in indices_with_data[1:]:
            #curr_flux_36 = self.band_36[index].flux_Jy()
            #curr_flux_45 = self.band_45[index].flux_Jy()

            curr_flux_36 = self.band_36[index].aperture_flux
            curr_flux_45 = self.band_45[index].aperture_flux

            self.mass_profile.append(\
                (arcsec_to_distance_ratio(index),\
                 flux_to_mstar(curr_flux_36, curr_flux_45, self.distance)))

    def radial_range_mass(self, low_limit, high_limit):
        return self.cumulative_mass_function(high_limit) - self.cumulative_mass_function(low_limit)

    def count_total_sne(self):
        if self.supernovae == None:
            raise AttributeError('Cannot count supernovae: this galaxy has not yet been cross-matched')
        else:
            return len(self.supernovae)

    def count_outskirts_sne(self):
        if self.supernovae == None:
            raise AttributeError('Cannot count supernovae: this galaxy has not yet been cross-matched')
        else:
            toReturn = 0
            for sn in self.supernovae:
                if sn.distance_ratio('read') > 1.0:
                    toReturn = toReturn + 1
            return toReturn

    # returns the average of the control times for total supernova rate calculation
    def mean_control_time(self):
        return (self.tc_Ia + self.tc_SE + self.tc_II) / 3.0

    def has_sfr(self):
        return (self.sfr_galspec != None) or (self.sfr_petrosian != None) or (self.sfr_SDSS != None) or (self.sfr_sersic != None)

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