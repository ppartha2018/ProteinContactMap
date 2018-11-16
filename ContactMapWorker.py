import sys
import os
import time
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
import numpy as np
#from display_contact_map import prepare_for_display
from compare_map import compare_contact_map
from Extract_Contact_Map import extract_raptorX_contact_map
from Extract_Contact_Map import construct_predicted_contact_map
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from ContactMapMetrics import ContactMapMetrics

warnings.simplefilter('ignore', PDBConstructionWarning)

#the Epsilon parameter - radius to be considered for nearest neighbor search
SEARCH_RADIUS = 8.0
# atoms to consider (CA or CB) for nearest neighbor calculations - eg: RaptorX is based in C-Beta, hence we compute for CB
ATOM_TO_CONSIDER = "CB" 
protein_id = ""
#output_prefix_path = "/scratch/ppartha/cm_result/"
output_prefix_path = "./cm_result_native_101/"
#predicted_cm_path_prefix = "/home/ppartha/shehu_lab/contact_map/RaptorX/"
predicted_cm_path_prefix = "./RaptorX/"


class ContactMapWorker:

    protein_id = ""
    pdb_location = ""
    contactMapMetrics = None

    def __init__(self, protein_id, pdb_location):
        self.protein_id = protein_id
        self.pdb_location = pdb_location
        #singleton object - one instance for a single pdb
        self.contactMapMetrics = ContactMapMetrics()
    
    '''
    Creates and returns the contact_map of atoms.
    The input is list of atom, usually c-alpha or c-beta atoms in the amino acid chain.
    The contact_map is a matrix MxM where M is size of atom list. 
    Each entry of contact_map C[i][j] = 1, if R[i][j] < Rc, 0 otherwise. Rc -> minimum threshold, distance between atoms given in angstrom unit.
    '''
    def get_original_contacts(self,atoms = [], residue_length = 1):
        neighbor_search = NeighborSearch(atoms)
        neighbors = neighbor_search.search_all(radius=SEARCH_RADIUS)
        #print neighbors
        #print("Size of Neighbhor Map: %d" % len(neighbors))
        #neighbors_map: dictionary with tuple(ca_i,ca_j) as keys - eg: [(ca1,ca2):1,(ca21,ca45):1,...]
        #In following line of code: x[0].get_parent() => atom.get_parent() which is residue. then, residue.get_id() is a tuple like ("", 10, ""),
        # continuing... here 10 is the residue's sequence identifier and thats what we want to use to define residue residue contacts..
        neighbors_map = {(x[0].get_parent().get_id()[1],x[1].get_parent().get_id()[1]): 1 for x in neighbors}
        no_of_original_contacts = len(neighbors_map) 
        #print neighbors_map
        print("In get_original_contacts: total no of residues: %d, no of original contacts: %d" % (residue_length,len(neighbors_map)))
        contact_map = {}
        #constuct contact map
        #total n(n-1)/2 entries, ( which excludes diagonal and symmetric entries from n*n pairs)
        # in range function start is inclusive, end is exclusive, so range(0, n-1) gives values from 0,1,2.. upto n-2
        for x in range(1, residue_length):
            for y in range(x+1, residue_length+1):
                #put {(Ci,Cj):1} if present in neighbors map, 0 otherwise
                #atoms[x].get_parent().get_id()[1] = residue id where the atom belongs to
                #** Wront implementation line below - don't need to go through atom's residue way, we already have residues based map
                #contact_map[(x,y)] = neighbors_map.get(((atoms[x].get_parent().get_id()[1],atoms[y].get_parent().get_id()[1])),0)
                contact_map[(x,y)] = neighbors_map.get((x,y),0)
        print("In get_original_contacts, contact map constructed:: Length of the contact map: %d" % (len(contact_map)))
        #print contact_map
        return contact_map, no_of_original_contacts  

    def process_pdb(self):
        '''
        since the input pdb file may be huge with thousands of models,
        read one model at a time and give to the function that uses biopython. 
        Biopython apparently tries to load entire pdb file into memory before processing thinking that 
        that pdb file contains only one model, we don't want that.
        '''
        #pdb_file = "pdb/"+protein_id+".pdb"
        output_path = output_prefix_path +"/" + self.protein_id + "/"
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        predicted_cm_path = predicted_cm_path_prefix + self.protein_id + "/contactmap.txt"
        print("Working with PDB: %s" % (self.pdb_location))
        
        model_content = []
        with open(self.pdb_location) as ip:
            for line in ip:
                model_content.append(line)
                if(line[0] == 'T'):
                    self.extract_atoms(protein_id, model_content, predicted_cm_path, output_path)
                    #clear the list
                    del model_content[:]
        ip.close()

    '''
    Read the PDB and update contact map for each model in a strcuture.
    Uses BioPython's methods to read and manipulate pdb.
    pdb file name is given as command-line argument.
    '''     
    def extract_atoms(self,protein_id, model_content, predicted_cm_path, output_path):
        #print "CAlling..."
        temp_pdb_file = output_path+"_temp_pdb.pdb"
        
        str1 = ""
        with open(temp_pdb_file,"w") as testing:
            testing.write(str1.join(model_content))
        testing.close()
        parser = PDBParser(PERMISSIVE=1)
        structure  = parser.get_structure("st",temp_pdb_file)
        #to store contact_map for each model
        contact_map_holder = {}
        residue_length = 0
        for model in structure.get_list():
            c_backbone_atoms = []
            #print "Total no of chains: ", len(model.get_list())
            for chain in model.get_list():
                residue_length = len(chain.get_list())
                for residue in chain.get_list():
                    #because GLY does not have CB atoms
                    if residue.get_resname() == "GLY":
                        c_backbone_atoms.append(residue["CA"])
                    elif residue.has_id(ATOM_TO_CONSIDER):
                        #can improve - need only the coordinates x,y,z - or not??
                        #print residue["CA"].get_serial_number()
                        c_backbone_atoms.append(residue[ATOM_TO_CONSIDER])
                        #print(c_backbone_atoms.get_coord())
        cm_result = self.get_original_contacts(c_backbone_atoms, residue_length) 
        #stores contact map for each model
        contact_map_holder[model.get_id()] = cm_result[0]
        no_of_original_contacts = cm_result[1]
        #prepare_for_display(contact_map[model.get_id()],model.get_id(), no_of_atoms)
        #display_contact_map(model.get_id(),contact_map[model.get_id()])
        #contact_comparison_scores = []
        predicted_contacts = extract_raptorX_contact_map(predicted_cm_path, residue_length)
        #top10
        print "Length of computed contact map: ", len(contact_map_holder[model.get_id()])
        print "Length of predicted contacts: ", len(predicted_contacts)
        print "Top10:"
        predicted_contact_map = construct_predicted_contact_map(predicted_contacts, 10, residue_length)
        contact_comparison_scores = self.contactMapMetrics.compare_contact_map(contact_map_holder[model.get_id()], predicted_contact_map, residue_length, no_of_original_contacts)
        precision_file = output_path+self.protein_id+"_cm_precision_top10.txt"
        coverage_file = output_path+self.protein_id+"_cm_coverage_top10.txt"
        recall_file = output_path+self.protein_id+"_cm_recall_top10.txt"
        f1_file = output_path+self.protein_id+"_cm_f1_top10.txt"

        self.write_results_to_file(precision_file, str(contact_comparison_scores[0]))
        self.write_results_to_file(coverage_file, str(contact_comparison_scores[1]))
        self.write_results_to_file(recall_file, str(contact_comparison_scores[2]))
        self.write_results_to_file(f1_file, str(contact_comparison_scores[3]))
        #top L/5
        print "TopL/5:"
        predicted_contact_map = construct_predicted_contact_map(predicted_contacts, residue_length/5, residue_length)
        contact_comparison_scores = self.contactMapMetrics.compare_contact_map(contact_map_holder[model.get_id()], predicted_contact_map, residue_length, no_of_original_contacts)
        precision_file = output_path+self.protein_id+"_cm_precision_topLBy5.txt"
        coverage_file = output_path+self.protein_id+"_cm_coverage_topLBy5.txt"
        recall_file = output_path+self.protein_id+"_cm_recall_topLBy5.txt"
        f1_file = output_path+self.protein_id+"_cm_f1_topLBy5.txt"

        self.write_results_to_file(precision_file, str(contact_comparison_scores[0]))
        self.write_results_to_file(coverage_file, str(contact_comparison_scores[1]))
        self.write_results_to_file(recall_file, str(contact_comparison_scores[2]))
        self.write_results_to_file(f1_file, str(contact_comparison_scores[3]))

        #top L/2
        print "TopL/2:"
        predicted_contact_map = construct_predicted_contact_map(predicted_contacts, residue_length/2, residue_length)
        contact_comparison_scores = self.contactMapMetrics.compare_contact_map(contact_map_holder[model.get_id()], predicted_contact_map, residue_length, no_of_original_contacts)
        precision_file = output_path+self.protein_id+"_cm_precision_topLBy2.txt"
        coverage_file = output_path+self.protein_id+"_cm_coverage_topLBy2.txt"
        recall_file = output_path+self.protein_id+"_cm_recall_topLBy2.txt"
        f1_file = output_path+self.protein_id+"_cm_f1_topLBy2.txt"

        self.write_results_to_file(precision_file, str(contact_comparison_scores[0]))
        self.write_results_to_file(coverage_file, str(contact_comparison_scores[1]))
        self.write_results_to_file(recall_file, str(contact_comparison_scores[2]))
        self.write_results_to_file(f1_file, str(contact_comparison_scores[3]))

    def write_results_to_file(self,filename, value):
        with open(filename, "a") as res_file:
            res_file.write(value)
            res_file.write("\n")
        res_file.close()
    
    def generate_contact_map_metrics(self):
        self.process_pdb()
    
    def get_precision(self):
        self.contactMapMetrics.get_precision()
        print("calculated precision:", self.contactMapMetrics.get_precision())
    
    def get_coverage(self):
        self.contactMapMetrics.get_coverage()
        print("calculated coverage:", self.contactMapMetrics.get_coverage())

# Program Arguments: Arg1: protein_id, Arg2: pdb location
# Usage: python contact_map2.py 1DTJA /home/xxx/yyy/pdb/1DTJA.pdb   
if __name__ == "__main__":
    #main(sys.argv[1])
    if(len(sys.argv)< 2):
	    print "Missing Arguments. Usage: python contact_map2.py 1AIL pdb_location"
    else:
        cmWorker = ContactMapWorker(sys.argv[1], sys.argv[2])
    	cmWorker.process_pdb()
        #cmWorker.get_precision()
        #cmWorker.get_coverage()