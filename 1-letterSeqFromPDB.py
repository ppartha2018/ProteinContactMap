"""
Created on Tue Jul 24 00:15:33 2018
@author: prasanna
"""
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import CaPPBuilder
'''
Extract the 1-letter amino acid sequence form the model. Only the first model is considered for
extracting the sequence since all of the models have the same sequence, thus efficient.
Requires BioPython Library 1.72.
'''
def extract_seq_from_models(protien_name, file_name, fasta_to_write):
    print "Working with.....: ", file_name
    #GET the first model from the pdb and write to a temp file
    #the code based on BioPython lib to extract the 1-letter sequence works fast if we just one model from a large pdb
    #thus reducing the memory
    temp_pdb_file = "_temp_pdb.pdb"
    with open(temp_pdb_file, "a") as temp:
        with open(file_name) as ip:
            for line in ip:
                temp.write(line)
                if(line[0] == 'T'):
                    break;
                
    structure = PDBParser().get_structure(protien_name,temp_pdb_file)
    # Using CA-CA
    ppb=CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        seq =  pp.get_sequence()

    os.remove(fasta_to_write)
    with open(fasta_to_write, "a") as fasta:
        fasta.write(">"+protien_name+"_seq\n")
        fasta.write(str(seq))
    ip.close()
    temp.close()
    fasta.close()
    os.remove(temp_pdb_file)
    print "Done, output file stored at: ", fasta_to_write
if __name__ == "__main__":
    extract_seq_from_models("1DTJA", "pdb/1DTJA_Output_RosettaDecoys_SubDirs0-99999.pdb", "fasta/1DTJA_Output_RosettaDecoys_SubDirs0-99999.fasta");