import sys
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
import numpy as np
#from display_contact_map import prepare_for_display
from compare_map import compare_contact_map
from Extract_Contact_Map import extract_raptorX_contact_map
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)

#the Epsilon parameter - radius to be considered for nearest neighbor search
SEARCH_RADIUS = 8.0
protein_id = ""
output_prefix_path = "/scratch/ppartha/"

'''
Creates and returns the contact_map of atoms.
The input is list of atom, usually c-alpha or c-beta atoms in the amino acid chain.
The contact_map is a matrix MxM where M is size of atom list. 
Each entry of contact_map C[i][j] = 1, if R[i][j] < Rc, 0 otherwise. Rc -> minimum threshold, distance between atoms given in angstrom unit.
'''
def get_contacts(atoms = []):
    neighbor_search = NeighborSearch(atoms)
    neighbors = neighbor_search.search_all(radius=SEARCH_RADIUS)
    #print("Size of Neighbhor Map: %d" % len(neighbors))
    #neighbors_map: dictionary with tuple(ca_i,ca_j) as keys - eg: [(ca1,ca2):1,(ca21,ca45):1,...]
    neighbors_map = {(x[0].get_serial_number(),x[1].get_serial_number()): 1 for x in neighbors} 

    contact_map = {}
    n = len(atoms)
    #constuct contact map
    #total n(n-1)/2 entries, ( which excludes diagonal and symmetric entries from n*n pairs)
    for x in range(1, n-1):
        for y in range((x+1),n):
            #put {(Ci,Cj):1} if present if present in neighbors map, 0 otherwise
            contact_map[(x,y)] = neighbors_map.get((atoms[x].get_serial_number(),atoms[y].get_serial_number()),0)
   # print("In get_contacts:Size of atoms: %d :Length of the contact map: %d" % (n,len(contact_map)))
    #print contact_map
    return contact_map  

def main(protein_id, predicted_cm_path):
    '''
    since the input pdb file may be huge with thousands of models,
    read one model at a time and give to the funciton that uses biopython. 
    Biopython apparently tries to load entire pdb file into memory before processing thinking that 
    that pdb file contains only one m
    odel, we don't want that.
    '''
    pdb_file = "pdb/"+protein_id+".pdb"
    print("Working with PDB: %s" % (pdb_file))
    
    model_content = []
    with open(pdb_file) as ip:
        for line in ip:
            model_content.append(line)
            if(line[0] == 'T'):
                extract_atoms(protein_id, model_content, predicted_cm_path)
                #clear the list
                del model_content[:]
    ip.close()
'''
Read the PDB and update contact map for each model in a strcuture.
Uses BioPython's methods to read and manipulate pdb.
pdb file name is given as command-line argument.
'''     
def extract_atoms(protein_id, model_content, predicted_cm_path):
    #print "CAlling..."
    temp_pdb_file = output_prefix_path+"_temp_pdb.pdb"
    cm_score1_file = output_prefix_path+protein_id+"_cm_score1.txt"
    cm_score2_file = output_prefix_path+protein_id+"_cm_score2.txt"
    cm_score3_file = output_prefix_path+protein_id+"_cm_accuracy.txt"
    
    str1 = ""
    with open(temp_pdb_file,"w") as testing:
        testing.write(str1.join(model_content))
    testing.close()
    parser = PDBParser(PERMISSIVE=1)
    structure  = parser.get_structure("st",temp_pdb_file)
    #to store contact_map for each model
    contact_map = {}
    for model in structure.get_list():
        c_alpha_atoms = []
        #print "Total no of chains: ", len(model.get_list())
        for chain in model.get_list():
            #print "Total no of residues: ", len(chain.get_list())
            for residue in chain.get_list():
                if residue.has_id("CA"):
                    #can improve - need only the coordinates x,y,z - or not??
                    #print residue["CA"].get_serial_number()
                    c_alpha_atoms.append(residue["CA"])
                    #print(c_alpha_atoms.get_coord())
	contact_map[model.get_id()] = get_contacts(c_alpha_atoms)
	#print("********Length of the contact map for model", model.get_id(), " is ", len(contact_map[model.get_id()]), "size of the residue chain: ", len(c_alpha_atoms))
        no_of_atoms = len(c_alpha_atoms)
    	#prepare_for_display(contact_map[model.get_id()],model.get_id(), no_of_atoms)
    	#display_contact_map(model.get_id(),contact_map[model.get_id()])
        scores = compare_contact_map(contact_map[model.get_id()],extract_raptorX_contact_map(predicted_cm_path,no_of_atoms),no_of_atoms)
        write_results_to_file(cm_score1_file, str(scores[0]))
        write_results_to_file(cm_score2_file, str(scores[1]))
        write_results_to_file(cm_score3_file, str(scores[2]))
        print "To write: %f in %s" % (scores[2],cm_score3_file)
        
def write_results_to_file(filename, value):
    with open(filename, "a") as res_file:
        res_file.write(value)
        res_file.write("\n")
	res_file.close()
        
        
   
if __name__ == "__main__":
    #main(sys.argv[1])
    if(len(sys.argv)< 2):
	print "Minium 1 argument needed. Usage: python contact_map2.py 1DTJ"
    else:
    	main(sys.argv[1], sys.argv[2])