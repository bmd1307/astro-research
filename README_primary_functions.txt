*************************************
*** Primary Function Descriptions ***
*************************************

********************************************************************************

plot_sne_cumulative()
	
	Description: Recreates the cumulative distribution function as seen in 
	Figure 8 of the Leaman paper (2011)

	Parameters: None

********************************************************************************

calc_sne_rate(n_bins)
	
	Description: 
	Recreates the (overall) supernova rate calculations as seen in Leaman and 
	Graur.

	Parameters:
	n_bins - the number of bins to calculate the rate

********************************************************************************

hist_total_sne_by_stellar_mass(n_bins = 10)
	
	Description: 
	Plots the step function (histogram) of the number of supernovae(total and 
	outskirts) by stellar mass.

	Parameters:
	n_bins - the number of bins for the step function (default 10)

********************************************************************************

total_sn_rate_outskirts(n_bins, 
						save_graph = False, 
						verbose = True, 
						show_graph = True,
						rate_function = sn_rate_outskirts, 
						title = 'Supernova Rate in Outskirts',
						yrange = [0.001, 30])
	
	Description: 
	Calculates the supernova rate in the outskirts. Creates plots similar to the
	rate plots in Graur's paper.The rate is calculated with fixed bins and a 
	sliding bins. The number of fixed bins is n_bins, and the sliding bin is 
	calculated with the same width of the fixed bins.

	Parameters:
	n_bins - the number of bins to calculate the rate
	save_graph - if True, saves the graph directly to a png file. File name = 
	out_rate_bins<n_bins>.png. False by default
	verbose - prints a lot of extra information about how the calculation is 
		progressing and detailed results. True by default
	show_graph - displays the graph and pauses the program. (default is True)
	rate_function - the function which performs the rate calculation valid 
		functions: sn_rate_outskirts (default), sn_rate_total
	title - the title for the graph. '(bins = <n_bins>)' is added to this name.
		default: 'Supernova Rate in Outskirts'
	yrange - the lower and upper limits of the y axis (default: [0.001, 30])

********************************************************************************

total_rate_dwarfs()

	Description: Calculates the supernova rate (number per millenium) in dwarf 
	galaxies (stellar mass < 109 solar masses). Details of the calculation are 
	printed.

	Parameters: None

********************************************************************************

total_rate_outskirts_spirals()

	Description: Calculates the supernova rate (number per millenium) in the 
	outskirts of spiral galaxies (r < R25). Details of the calculation are 
	printed.

	Parameters: None

********************************************************************************

compare_outskirts_to_dwarfs()

	Description: Runs total_rate_dwarfs() and total_rate_ouskirts_spirals() for 
	comparing these rates.

	Parameters: None

********************************************************************************

sne_radial_data()

	Description: Finds the radial dependence of supernovae. Groups the 
	supernovae into an ideal number of bins by radius, and fits a power function
	to the supernova density as a function of radius. Displays a graph of the 
	supernova density as a function of galactocentric radius.

	Parameters: None

********************************************************************************

sne_radial_histogram()

	Description: Plots a histogram of the number of supernovae vs radius. The 
	number of bins is chosen by the Freedman-Diaconis equation.

	Parameters: None

********************************************************************************

__main__()

	The First function called. Modify this to execute the desired function.
