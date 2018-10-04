import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.stats.stats import pearsonr
from matplotlib.patches import Rectangle

def read_from_file(file_name):
    values = []
    with open(file_name) as fp:
            for line in fp:
                values.append(float(line))
    return values   

def do_scatter_plot(measure1, measure2, measure1_file, measure2_file, path_prefix, protein_id):
    #test data
    #rmsd = [2.9680116920666535, 3.2662908745233143, 1.2865403599752712, 10.469770730882399, 5.514320965750489]
    #cm_accuracy = [0.277108433735, 0.316417910448, 0.251552795031, 0.311178247734, 0.339031339031]

    try:
        print("plotting for protien: %s measure1: %s, measure2: %s" % (pid, measure1, measure2))
        file_name = path_prefix + protein_id + "_" + measure1 + "_" + measure2 + ".png"
        m1 = read_from_file(measure1_file)
        m2 = read_from_file(measure2_file)

        plt.xlabel(measure1)
        plt.ylabel(measure2)
        plt.scatter(m1, m2)
        pearson_correlation = (pearsonr(m1, m2))[0]
        text1 = ("%s - Pearson Correlation: %s" % (pid, str(round(pearson_correlation,2))))
        extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        plt.legend([extra],[text1], loc='upper right')
        plt.savefig(file_name)
        plt.gcf().clear()
    except:
        print("Unexpected error:", sys.exc_info()[0])
        plt.gcf().clear()
   
#do_scatter_plot("contact_accuracy","rmsd","cm_scores/accuracy/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_accuracy.txt","rmsd/1DTJA_rmsd.txt")
#do_scatter_plot("contact_accuracy","score12","cm_scores/accuracy/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_accuracy.txt","energy_score12/1DTJA_Output_RosettaDecoys_SubDirs0-99999.score12")
#do_scatter_plot("contact_coverage","rmsd","cm_scores/coverage/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_coverage.txt","rmsd/1DTJA_rmsd.txt", "plots/", "1DTJA")
#do_scatter_plot("contact_coverage","score12","cm_scores/coverage/1DTJA_Output_RosettaDecoys_SubDirs0-99999_cm_coverage.txt","energy_score12/1DTJA_Output_RosettaDecoys_SubDirs0-99999.score12", "plots/", "1DTJA")

output_path_prefix = "cm_result_0929/plots/"

#os.makedirs(output_path_prefix)
#complete list - not using now as we don,t have rmsds for everything
#protein_ids = ["1AIL", "1AOY", "1C8CA", "1CC5", "1DTJA", "1DTDB", "1FWP", "1HHP", "1HZ6A", "1ISUA", "1SAP", "1TIG", "1WAPA", "2CI2", "2EZK", "2H5ND", "2HG6"]
protein_ids = ["1AIL", "1AOY", "1C8CA", "1DTJA", "1DTDB", "1HZ6A", "1ISUA", "1SAP", "1TIG", "2EZK"]
#protein_ids = ["1AIL", "1DTDB"]
for pid in protein_ids:
    precision_file = "cm_result_0929/" + pid + "/" + pid + "_cm_precision_topLBy5.txt"
    #coverage_file = "cm_result_0929/" + pid + "/" + pid + "_cm_coverage.txt"
    rmsd_file = "rmsd/" + pid + "_rmsd.txt"
    energy_file = "energy/" + pid + ".score12" 
    do_scatter_plot("contact_precision","rmsd", precision_file, rmsd_file, output_path_prefix, pid)
    do_scatter_plot("contact_precision","score12", precision_file, energy_file, output_path_prefix, pid)
    #do_scatter_plot("contact_coverage","rmsd", coverage_file, rmsd_file, output_path_prefix, pid)
    #do_scatter_plot("contact_coverage","score12", coverage_file, energy_file, output_path_prefix, pid)