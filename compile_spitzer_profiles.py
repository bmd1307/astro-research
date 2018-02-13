# compiles the radial spitzer profiles

import os
import csv

print(os.getcwd())

out_file = open('radial_profiles.csv', 'w', newline = '')

out_writer = csv.writer(out_file)

ellipse_fits_file_names = os.listdir('ellipse_fits')

out_writer.writerow(['Galaxy name',\
                    'Band frequency (10^-1 um)',\
                    'Ellipticity',\
                    'Ellipticity error',\
                    'Equivalent radius (arcsec)',\
                    'Flux in isophote (MJy/sr)',\
                    'Semimajor axis radius (arcsec)'])

counter = 1

for fname in ellipse_fits_file_names:

    curr_file = open('ellipse_fits\\' + fname)

    for line in curr_file:

        line_vals = line.split()

        name_params = fname.split('_')

        out_writer.writerow([name_params[0],\
                              name_params[1][:2],\
                              line_vals[6 - 1],\
                              line_vals[7 - 1],\
                              line_vals[41 - 1],\
                              line_vals[44 - 1],\
                              line_vals[51 - 1]])
    if counter % 10 == 0:
        print('.', end = '')
        if counter % 800 == 0:
            print()
    counter = counter + 1

out_file.close()
