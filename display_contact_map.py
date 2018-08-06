import numpy as np
import matplotlib.pyplot as plt
import time

def plot_matrix(data):
    fig = plt.figure(dpi=150)
    plt.imshow(data)
    plt.colorbar()
    plt.show()
    fig.savefig("contact_map.png");

#contact_map is of type {(Ci,Cj):1/0,...}, i,j ranges from 1 to n
#convert it into a symmetric matrix and feed it into the display method
def prepare_for_display(contact_map,model_id,length_of_residue_chain):
    # 0 to n-1 rows, 0 to n-1 columns
    contact_map_matrix = np.zeros((length_of_residue_chain,length_of_residue_chain))
    for x,y in contact_map.iteritems():
        if(y == 1):
            contact_map_matrix[x[0]][x[1]] = 1
            contact_map_matrix[x[1]][x[0]] = 1
    #display_contact_map(model_id,contact_map_matrix)
    plot_matrix(contact_map_matrix)


