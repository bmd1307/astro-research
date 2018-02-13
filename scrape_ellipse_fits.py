
from bs4 import BeautifulSoup
import urllib.request


irsa_url = 'http://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/'

dest_dir = 'C:\\Users\\Brennan\\Desktop\\astro-research\\ellipse_fits'

galaxy_list_page = urllib.request.urlopen(irsa_url)

s = str(galaxy_list_page.read())

b = BeautifulSoup(s, 'html.parser')

links = [link.string[:-1] for link in b.body.pre.find_all('a') if link.string[-1] == '/']

counter = 1

print('Parsing')
for l in links:
    print('.', end = '')
    if counter % 80 == 0:
        print()
    file_3_6_band = irsa_url + l + '/P3/' + l + '.1fr6a_noclean_fin.dat'
    file_4_5_band = irsa_url + l + '/P3/' + l + '.2fr6a_noclean_fin.dat'

    fname_36 = l + '_' + '36u.dat'
    fname_45 = l + '_' + '45u.dat' 

    try:
        urllib.request.urlretrieve(file_3_6_band, dest_dir + '\\' + fname_36)
        urllib.request.urlretrieve(file_4_5_band, dest_dir + '\\' + fname_45)
    except urllib.error.HTTPError as e:
        if e.getcode() == 404:
            continue
        else:
            raise

    
print('\nDone')
print('Parsed ellipse data for', counter, 'galaxies.')
