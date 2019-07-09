def crossing_finder(x_lst, y_lst):
    master_x_cross = []
    master_y_cross = []
    for i in range(len(y_lst)):
        if y_lst[i]*y_lst[i+1] <= 0:
            x_cross_lst = [x_lst[i-1], x_lst[i], x_lst[i+1], x_lst[i+2]]  
            y_cross_lst = [y_lst[i-1], y_lst[i], y_lst[i+1], y_lst[i+2]]
            master_x_cross.append(x_cross_lst)
            master_y_cross.append(y_cross_lst)
    return master_x_cross, master_y_cross


