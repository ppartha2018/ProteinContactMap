import numpy as np
import matplotlib.pyplot as plt
import os

def read_from_file(file_name):
    values = []
    with open(file_name) as fp:
            for line in fp:
                values.append(float(line))
    return values

def do_scatter_plot(measure1, measure2, measure1_file, measure2_file):
    #test data
    #rmsd = [2.9680116920666535, 3.2662908745233143, 1.2865403599752712, 10.469770730882399, 5.514320965750489]
    #cm_accuracy = [0.277108433735, 0.316417910448, 0.251552795031, 0.311178247734, 0.339031339031]

    m1 = read_from_file(measure1_file)
    m2 = read_from_file(measure2_file)

    plt.xlabel(measure1)
    plt.ylabel(measure2)
    plt.scatter(m1, m2)
    plt.show()

#do_scatter_plot("contact_accuracy","rmsd","cm_scores/accuracy/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_accuracy.txt","rmsd/1DTJA_rmsd.txt")
#do_scatter_plot("contact_accuracy","score12","cm_scores/accuracy/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_accuracy.txt","energy_score12/1DTJA_Output_RosettaDecoys_SubDirs0-99999.score12")
do_scatter_plot("contact_coverage","rmsd","cm_scores/coverage/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_coverage.txt","rmsd/1DTJA_rmsd.txt")
do_scatter_plot("contact_coverage","score12","cm_scores/coverage/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_coverage.txt","energy_score12/1DTJA_Output_RosettaDecoys_SubDirs0-99999.score12")