

sne_file = open('sn-full.txt')

sne_types = []

for line in sne_file:
    curr_val = line[87:87+7]

    sne_types.append(curr_val.strip())

sne_file.close()

sne_types = sorted(sne_types)

type_hist = {}

for snt in sne_types:
    if snt in type_hist:
        type_hist[snt] = type_hist[snt] + 1
    else:
        type_hist[snt] = 1

sne_keys = type_hist.keys()

sne_keys = sorted(sne_keys)

type_II_keys = [k for k in sne_keys if 'II' == k[0:2]]

type_Ia_keys = [k for k in sne_keys if 'Ia' == k[0:2]]

type_Ib_keys = [k for k in sne_keys if 'Ib' == k[0:2] and 'Ibc' != k[0:3]]

type_Ibc_keys = [k for k in sne_keys if 'Ibc' == k[0:3]]

type_Ic_keys = [k for k in sne_keys if 'Ic' == k[0:2]]

print('*** FORMAT ***')
print('General SN Type: <# General SNe>')
print('\tSpecific SN Type: <# Specific SNe>')

print('*** DATA ***')

print('Type Ia: ', sum([type_hist[k] for k in type_Ia_keys]))
for k in type_Ia_keys:
    print('\t', (k + ':').ljust(8), str(type_hist[k]).rjust(4), sep='')

print('Type Ib: ', sum([type_hist[k] for k in type_Ib_keys]))
for k in type_Ib_keys:
    print('\t', (k + ':').ljust(8), str(type_hist[k]).rjust(4), sep='')

print('Type Ibc: ', sum([type_hist[k] for k in type_Ibc_keys]))
for k in type_Ibc_keys:
    print('\t', (k + ':').ljust(8), str(type_hist[k]).rjust(4), sep='')

print('Type Ic: ', sum([type_hist[k] for k in type_Ic_keys]))
for k in type_Ic_keys:
    print('\t', (k + ':').ljust(8), str(type_hist[k]).rjust(4), sep='')

print('Type II: ', sum([type_hist[k] for k in type_II_keys]))
for k in type_II_keys:
    print('\t', (k + ':').ljust(8), str(type_hist[k]).rjust(4), sep='')

