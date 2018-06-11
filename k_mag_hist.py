# k mag histogram

import matplotlib.pyplot as plt
import math
import scipy
import scipy.integrate

from scipy.special import gammaincc

M_star = -20.83
L_star = 2e10

def lum_density(La, Lb):
    phi_star = 1.46e-2
    L_star = 2e10
    alpha = -1.2

    return phi_star * L_star * (gammaincc(alpha + 2, La / L_star) - gammaincc(alpha + 2, Lb / L_star))

# calculates phi(lum) for a luminosity in solar masses
def phi(lum):
    phi_star = 1.46e-2
    alpha = -1.2
    return phi_star * (lum/L_star)**alpha * math.exp(-(lum/L_star)) / L_star

def schechter(low_lum, high_lum):
    #return coef * (gammaincc(alpha + 1, low_lum) - gammaincc(alpha + 1, high_lum))
    return scipy.integrate.quad(phi, low_lum, high_lum)

def schechter_density(low_lum, high_lum):
    return scipy.integrate.quad((lambda l: (l) * phi(l)), low_lum, high_lum)

def schechter_mag(low, high):
    coef = 1.46e-2
    alpha = -1.20

    lum_f = (lambda m: 10**(0.4 * (M_star - m)))
    schechter_f = (lambda m: 0.4 * math.log10(10) * coef * \
                   lum_f(m) ** (alpha + 1) * \
                   math.exp(-lum_f(m)))
    return scipy.integrate.quad(schechter_f, low, high)
    

gal_file = open('galaxy-full.txt')

graur_file = open('table2.dat')

k_mags = []
k_errs = []

k_lums = []
distances = []

graur_lums = []

for line in gal_file:
    vals = line[163:175].split()
    dist = float(line[67:75])

    if vals[0] == '99.999':
        continue
    
    k_mags.append(float(vals[0]))
    k_errs.append(float(vals[1]))
    distances.append(dist)

for line in graur_file:
    graur_lums.append(float(line[67:75] )* 10**10)

graur_lums = [lum for lum in graur_lums if lum != 0.0]

abs_mags = []

for i in range(0, len(distances)):
    curr_dist = distances[i] * 1e6 # convert to pc
    curr_mag = k_mags[i]

    curr_abs_mag = curr_mag - 5 * math.log10(curr_dist / 10)

    abs_mags.append(curr_abs_mag)


for curr_mag in abs_mags:
    k_lums.append(10 ** (0.4 * (M_star - curr_mag)))

print(graur_lums[0:20])

print('min mag:', min(abs_mags))
print('max mag:', max(abs_mags))

print('min lum:', min(graur_lums))
print('max lum:', max(graur_lums))

print('Schechter min-max:', schechter(min(graur_lums), max(graur_lums)))

max_dist = max(distances)

print('max dist', max_dist)

max_vol = (4/3) * math.pi * max_dist**3

print('Est total galaxies lum:', max_vol * schechter(min(graur_lums), max(graur_lums))[0])



bins = [10**(6+n/3) for n in range(0, 20)]
print(bins)

graur_hist = []
for i in range(len(bins) - 1):
    curr_count = 0
    for lum in graur_lums:
        if bins[i] < lum <= bins[i + 1]:
            curr_count = curr_count + 1
    graur_hist.append(curr_count)

print(graur_hist)

total_dens = schechter(bins[0], bins[-1])[0]

volume = (4/3) * math.pi * (60.0) ** 3

print('Total density:', '%1.3e' % total_dens)

for i in range(0, len(bins) - 1):
    print('%1.3e' % bins[i], '%1.3e' % bins[i + 1],\
          #lum_density(bins[i], bins[i+1]) / (0.5*bins[i] + 0.5*bins[i+1]),\
          schechter(bins[i], bins[i+1])[0],\
          '%15.3f' % (schechter(bins[i], bins[i+1])[0] / total_dens),\
          '%15.3f' % (schechter(bins[i], bins[i+1])[0] * volume),\
          '%8i' % graur_hist[i],\
          sep='\t')

print('*** DENSITY ***')
for i in range(0, len(bins) - 1):
    print('%1.3e' % bins[i], '%1.3e' % bins[i + 1],\
          '%15.3e' % (schechter_density(bins[i], bins[i+1])[0]),\
          sep='\t')

print('******')

print('Total Luminosity of Graur galaxies (Lsun):', '%1.3e' % sum(graur_lums))
print('Volume of graur sample (Mpc):', volume)
print('Graur overall luminosity density:', '%1.3e' % (sum(graur_lums) / volume))
print('Schechter density 1e6 to 1e12', '%1.3e' % schechter_density(1e6, 1e12)[0])
print('Inferred SN count from schechter density', '%8i' % (schechter_density(1e6, 1e12)[0] * 1.3e-2))

gal_file.close()
graur_file.close()
