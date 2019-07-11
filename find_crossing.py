def crossing_finder(x_lst, y_lst):
    master_x_cross = []
    master_y_cross = []
    cross_indicies = [] # where the crossings occur
    i = 1
    full_orbit = False
    while i < len(y_lst) - 2:#Makes sure it doesn't do an out of bounds index
            if y_lst[i]*y_lst[i+1] <= 0:
                if full_orbit:
                    x_cross_lst = [x_lst[i-1], x_lst[i], x_lst[i+1], x_lst[i+2]]  
                    y_cross_lst = [y_lst[i-1], y_lst[i], y_lst[i+1], y_lst[i+2]]
                    master_x_cross.append(x_cross_lst)
                    master_y_cross.append(y_cross_lst)
                    cross_indicies.append(i - 1) # add index of where the set of four starts

                # if y_lst[i+1] == 0:#Makes sure it doesn't add repeating sequences
                #     i += 2  
                #     print('second:',i)

                full_orbit = not full_orbit #Make sure it skips every other value
            i += 1
    return cross_indicies, master_x_cross, master_y_cross










