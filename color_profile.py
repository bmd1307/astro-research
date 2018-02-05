
class ColorData:
    def __init__(self, file_line):
        self.galaxy_name    = file_line[0:10].strip() # get rid of extra spaces
        self.radius_sec     = float(file_line[11:14])

        self.mu_FUV         = float(file_line[21:27])
        self.mu_FUV_err     = float(file_line[28:36])
        self.cumu_FUV       = float(file_line[37:43])
        self.cumu_FUV_err   = float(file_line[44:51])

        self.mu_NUV         = float(file_line[52:58])
        self.mu_NUV_err     = float(file_line[59:67])
        self.cumu_NUV       = float(file_line[68:74])
        self.cumu_NUV_err   = float(file_line[75:80])

        self.mu_36          = float(file_line[101:107])
        self.mu_36_err      = float(file_line[115:124])

    # add functions which convert a magnitude to stellar mass/frequency

def parse_color_profile_file(color_profile_dict, file_name):
    print('Parsing file ' + file_name)

    # color_profile_dict is a dictionary from galaxy names to a dictionary
    # of radii to color_data

    color_file = open(file_name)

    for line in color_file:
        # get each color data object
        curr_color_data = ColorData(line)

        if not (curr_color_data.galaxy_name in color_profile_dict.keys()):
            color_profile_dict[curr_color_data.galaxy_name] = {}

        color_profile_dict[curr_color_data.galaxy_name]\
                          [curr_color_data.radius_sec] = curr_color_data

    # print some information on number of