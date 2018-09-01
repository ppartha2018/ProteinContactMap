# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 00:15:33 2018

@author: prasanna
"""

#from display_contact_map import prepare_for_display

def extract_raptorX_contact_map(file_name, length_of_residue_chain):
    contact_map = {}
    cm_file = open(file_name,"r")
    n = length_of_residue_chain
    
    for line in cm_file:
       if(line[0].isdigit()):
           a_line = line.split()
           contact_map[(int(a_line[0]),int(a_line[1]))] = 1
    #print("no of contacts in predicted contact map: %d" % (len(contact_map))) 
    for x in range(n-1):
        for y in range((x+1),n):
            #put other entries to 0
            contact_map[(x,y)] = contact_map.get((x,y),0)
            
   # prepare_for_display(contact_map,1,length_of_residue_chain)
    return contact_map
 
if __name__ == "__main__":
    extract_raptorX_contact_map("RaptorX_Prediction_1DTJ\contactmap.txt",266);