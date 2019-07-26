import sys
import numpy as np


def main(argv):
    filename = argv[0]
    overlay = False
    print(argv[1])
    if argv[1] == '-ov' or argv[1] == '--overlay':
        overlay = True
        orbits_to_get = argv[2:]
        print('Calculating overlay.')
    else:
        orbits_to_get = argv[1:]
    

    data = np.genfromtxt(filename, names=True)

    orbit_list = []
    orbit_start_index = 0
    i = 0
    orbit_start_value = data['star_angle'][i]
    for i in range(len(data['star_angle'])):
        if data['star_angle'][i] > orbit_start_value + 2 * np.pi:
            orbit_list.append((orbit_start_index, i))
            orbit_start_index = i
            orbit_start_value = data['star_angle'][i]


    header_text = ''
    for s in data.dtype.names:
        header_text += s + ' '

    for orbit_num_str in orbits_to_get:
        orbit_num_int = int(orbit_num_str)

        outfile_name = "" # name it differently
        if orbit_num_int < 0: # no negative indecies
            orbit_num_int = len(orbit_list) + orbit_num_int
            outfile_name = str(orbit_num_int) + filename
        else:
            outfile_name = orbit_num_str + filename

        print('Getting orbit number', orbit_num_int)

        # name it differently for an overlay
        if overlay is True:
            outfile_name = 'ov' + outfile_name
        else:
            outfile_name = 'o' + outfile_name

        print('Trying to write to', outfile_name, end='')

        begin, end = orbit_list[orbit_num_int]
        sub_data = data[begin:end + 1]
        if overlay is True and orbit_num_int is not 0: # shift it to overlay
            sub_data['star_angle'] = [x - (2 * np.pi * orbit_num_int) for x in sub_data['star_angle']]
        np.savetxt(outfile_name, sub_data, header=header_text)
        print('...succeeded')


if __name__ == "__main__":
    main(sys.argv[1:])
