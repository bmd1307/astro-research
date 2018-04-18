import numpy as np

num_trials = 100000

def calc_types(num_Ia, num_SE, num_II):
    spiral_Ia = np.random.poisson(num_Ia, num_trials)
    spiral_SE = np.random.poisson(num_SE, num_trials)
    spiral_II = np.random.poisson(num_II, num_trials)

    ratios = [(SE+II)/Ia for (Ia, SE, II) in zip(spiral_Ia, spiral_SE, spiral_II) if Ia != 0]

    one_sig = np.std(ratios)

    print("Outskirts Spirals:")
    print("Type Ia:", num_Ia, "Type SE:", num_SE, "Type II:", num_II)
    print("Mean Ratio: ", "%1.3f" % np.mean(ratios))
    print("Stdev Ratio:", "%1.3f" % np.std(ratios))
    print("14th percentile:", "%1.3f" % np.percentile(ratios, 16.0, interpolation = 'linear'))
    print("50th percentile:", "%1.3f" % np.percentile(ratios, 50.0, interpolation = 'linear'))
    print("86th percentile:", "%1.3f" % np.percentile(ratios, 84.0, interpolation = 'linear'))


print('Total Galaxy Sample')
calc_types(372, 144, 399)
print()

print('Total Galaxy Outskirts')
calc_types(64, 11, 41)
print()

print('Spiral Outskirts')
calc_types(43, 11, 32)
print()

print('Dwarfs')
calc_types(5, 5, 6)
print()
