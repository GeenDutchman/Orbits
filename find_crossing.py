def crossing_finder(x_lst, y_lst):
    master_x_cross = []
    master_y_cross = []
    i = 1
    while i < len(y_lst) - 2:#Makes sure it doesn't do an out of bounds index
            if y_lst[i]*y_lst[i+1] <= 0:
                x_cross_lst = [x_lst[i-1], x_lst[i], x_lst[i+1], x_lst[i+2]]  
                y_cross_lst = [y_lst[i-1], y_lst[i], y_lst[i+1], y_lst[i+2]]
                master_x_cross.append(x_cross_lst)
                master_y_cross.append(y_cross_lst)
        
                if y_lst[i+1] == 0:#Makes sure it doesn't add repeating sequences
                    i += 2  
                    print('second:',i)
            i += 1
    return master_x_cross, master_y_cross

#Goals:
#Make sure it skips every odd value, only considers even values








