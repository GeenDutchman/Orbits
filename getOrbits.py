import sys
import numpy as np


def main(argv):
    filename = argv[0]
    outfile = argv[1]
    orbits_to_get = argv[2:]

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

    for orbit_num in orbits_to_get:
        outfile_name = 'o' + str(orbit_num) + outfile
        print('Trying to write to', outfile_name, end='')
        begin, end = orbit_list[int(orbit_num)]
        sub_data = data[begin:end + 1]
        np.savetxt(outfile_name, sub_data, header=header_text)
        print('...succeeded')


if __name__ == "__main__":
    main(sys.argv[1:])
