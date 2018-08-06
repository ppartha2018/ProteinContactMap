from __future__ import division
from Bio.PDB import *

def cmap_key_equals(key1,key2):
    #Key1, Key2 are tuples where each element is of type Bio.PDB.Atom.Atom
    return ((key1[0].__eq__(key2[0])) and (key1[1].__eq__(key2[1])))

'''
Contact map evaluation measures:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3226919/
'''
def compare_contact_map(computed_map = {}, predicted_map = {}, no_of_atoms = 1):
    #example
    #computed_map = {('ca1','ca2'):0,('ca1','ca1'):1,('ca1','ca3'):1}
    #predicted_map = {('ca1','ca2'):1,('ca1','ca1'):1,('ca1','ca3'):1}
    
    #print("Length of computed map: %d, Length of predicted map: %d" % (len(computed_map), len(predicted_map)))
    no_of_contacts = (no_of_atoms * (no_of_atoms - 1)) / 2
    computed_map_size = len(computed_map)
    predicated_map_size = len(predicted_map)
    computed_map_1s = 0
    computed_map_0s = 0
    predicted_map_1s = 0
    predicted_map_0s = 0
    true_positives = 0
    false_positives = 0
    
    for key,value in computed_map.iteritems():
        #print key,value
        if value == 1:
            computed_map_1s = computed_map_1s + 1
            if(value == predicted_map[key]):
                true_positives = true_positives + 1
        else:
            computed_map_0s = computed_map_0s + 1
            if(value != predicted_map[key]):
                false_positives = false_positives + 1
        #if(value == predicted_map[key]):
        #    difference = difference + 1;
    
    for key,value in predicted_map.iteritems():
        #print key,value
        if value == 1:
            predicted_map_1s = predicted_map_1s + 1;
        else:
            predicted_map_0s = predicted_map_0s + 1;
    
    #print computed_map_0s, computed_map_1s, predicted_map_0s, predicted_map_1s, difference
    score1 = predicted_map_1s / no_of_atoms;
    score2 = (predicted_map_0s + predicted_map_1s) / no_of_atoms;
    
    scores = []
    scores.append(score1)
    scores.append(score2)
    #Accuracy Acc = TP/(TP + FP), 
    #where TP and FP are the numbers of correctly and incorrectly predicted contacts, respectively
    accuracy = true_positives / (true_positives + false_positives);
    scores.append(accuracy)
    
    #print ("Evaluation: Score1: %f, Score2: %f, Score3: %f" % (score1,score2,score3))
    #print ("Evaluation: Accuracy: %f" % (accuracy))
    return scores
    
    #cmap_key_equals(computed_map[('ca1','ca2')],computed_map[]); 
    #print computed_map[('ca1','ca2')]

if __name__ == "__main__":
    compare_contact_map()
