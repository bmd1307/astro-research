# k mag histogram

import matplotlib.pyplot as plt

gal_file = open('galaxy-full.txt')

graur_file = open('table2.dat')

k_mags = []
k_errs = []

k_lums = []

for line in gal_file:
    vals = line[163:175].split()

    if vals[0] == '99.999':
        continue
    
    k_mags.append(float(vals[0]))
    k_errs.append(float(vals[1]))

for line in graur_file:
    value = float(line[67:75])

    k_lums.append(value)

print(k_mags[0:20])
print(k_lums[0:20])

plt.hist(k_mags, bins=20)
plt.show()

plt.hist(k_lums)
plt.yscale('log')
plt.show()

gal_file.close()
graur_file.close()
