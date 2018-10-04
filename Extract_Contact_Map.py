# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 00:15:33 2018

@author: prasanna
"""

#from display_contact_map import prepare_for_display

def extract_raptorX_contact_map(file_name, length_of_residue_chain):
    print("Extracting predicted contact map from: %s" % (file_name))
    cm_file = open(file_name,"r")
    print("In extract_contact_map: length of the residue chain: %d" % (length_of_residue_chain))
    predicted_contacts = []
    for line in cm_file:
       if(line[0].isdigit()):
            a_line = line.split()
            predicted_contacts.append((int(a_line[0]),int(a_line[1])))
    print("After Extracting: no of predicted contacts from file: ", len(predicted_contacts))
    return predicted_contacts

def construct_predicted_contact_map(predicted_contacts, limit, length_of_residue_chain):
    contact_map = {}
    print "length of residue chain: ", length_of_residue_chain
    for contact in predicted_contacts[:limit]:
        contact_map[contact] = 1         
    print("In construct_predicted_contact_map: No of predicted contacts in predicted contact map: %d, limit: %d" % (len(contact_map), limit))

    # beacause the residue ranges start from 1 upto n
    for x in range(1, length_of_residue_chain):
        for y in range((x+1), length_of_residue_chain+1):
            #put other entries to 0
            contact_map[(x,y)] = contact_map.get((x,y),0)

    print("In construct_predicted_contact_map: No of entries in predicted contact map after assigning for non contacts: %d" % (len(contact_map)))      
    #prepare_for_display(contact_map,1,length_of_residue_chain)
    return contact_map
 
if __name__ == "__main__":
    extract_raptorX_contact_map("RaptorX_Prediction_1DTJ\contactmap.txt",266);