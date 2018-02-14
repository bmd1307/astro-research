
import csv
import math

class SpitzerRecord:
    def __init__(self, line_data):
        # unlike the other __init__ functions, this one comes from a CSV file, so line_data is an array
        self.galaxy_name = line_data[0]

        self.band = line_data[1]

        self.ellipticity =              float(line_data[2]) if line_data[2] != 'INDEF' else None

        self.ellipticity_err =          float(line_data[3]) if line_data[3] != 'INDEF' else None

        self.equivalent_radius_sec =    float(line_data[4]) if line_data[4] != 'INDEF' else None

        self.aperture_flux =            float(line_data[5]) if line_data[5] != 'INDEF' else None # units: MJy / sr

        self.maj_axis_radius_sec =      int(float(line_data[6]))

    def flux_Jy(self):
        if self.aperture_flux is None or self.ellipticity is None or self.maj_axis_radius_sec is None:
            return None
        # converts:    MJy/sr -> Jy/arcsec^2                        *** ellipse area *** ((1 - eta) * pi * a^2
        return (2.35045e-5 * self.aperture_flux) * math.pi * (1 - self.ellipticity) * self.maj_axis_radius_sec ** 2

def parse_spitzer_file(spitzer_dict, file_name):
    print('Parsing file ' + file_name)

    spitzer_file = open(file_name)
    csv_reader = csv.reader(spitzer_file)

    # skip the header line
    next(csv_reader)

    for line in csv_reader:
        curr_spitzer_data = SpitzerRecord(line)

        if not (curr_spitzer_data.galaxy_name in spitzer_dict):
            spitzer_dict[curr_spitzer_data.galaxy_name] = ({}, {})

        if curr_spitzer_data.band == '36':
            spitzer_dict[curr_spitzer_data.galaxy_name][0][curr_spitzer_data.maj_axis_radius_sec] = \
                curr_spitzer_data
        elif curr_spitzer_data.band == '45':
            spitzer_dict[curr_spitzer_data.galaxy_name][1][curr_spitzer_data.maj_axis_radius_sec] = \
                curr_spitzer_data
        else:
            raise AttributeError('Invalid band width found while parsing spitzer data.')