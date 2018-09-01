from __future__ import division
from Bio.PDB import *

class ContactMapMetrics:

    total_samples = 0
    accuracy_sum = 0
    coverage_sum = 0
    accuracy_average = 0
    coverage_average = 0

    def cmap_key_equals(self,key1,key2):
    #Key1, Key2 are tuples where each element is of type Bio.PDB.Atom.Atom
        return ((key1[0].__eq__(key2[0])) and (key1[1].__eq__(key2[1])))
    
    '''
    Contact map evaluation measures:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3226919/
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1404-z
    ''' 
    def compare_contact_map(self, computed_map = {}, predicted_map = {}, no_of_atoms = 1, no_of_native_contacts = 1):
        #example
        #computed_map = {('ca1','ca2'):0,('ca1','ca1'):1,('ca1','ca3'):1}
        #predicted_map = {('ca1','ca2'):1,('ca1','ca1'):1,('ca1','ca3'):1}
        
        #print("Length of computed map: %d, Length of predicted map: %d" % (len(computed_map), len(predicted_map)))

        accuracy = 0
        coverage = 0

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
                if(value == predicted_map[key]): #if predicted as contact
                    true_positives = true_positives + 1
                else:
                    false_positives = false_positives + 1
            else:
                computed_map_0s = computed_map_0s + 1
            #if(value == predicted_map[key]):
            #    difference = difference + 1;
        
        for key,value in predicted_map.iteritems():
            #print key,value
            if value == 1:
                predicted_map_1s = predicted_map_1s + 1;
            else:
                predicted_map_0s = predicted_map_0s + 1;
        
        #print computed_map_0s, computed_map_1s, predicted_map_0s, predicted_map_1s, difference
        
        #Accuracy Acc = TP/(TP + FP), 
        #where TP and FP are the numbers of correctly and incorrectly predicted contacts, respectively
        accuracy = true_positives / (true_positives + false_positives)
        #coverage = TP/ Nc, Nc - no of contacts in native model
        coverage = predicted_map_1s / no_of_native_contacts
        #score2 = (predicted_map_0s + predicted_map_1s) / no_of_atoms
        
        scores = []
        scores.append(accuracy)
        scores.append(coverage)
        #scores.append(score2)
        
        
        #print ("Evaluation: Score1: %f, Score2: %f, Score3: %f" % (score1,score2,score3))
        #print ("Evaluation: Accuracy: %f" % (accuracy))
        self.total_samples = self.total_samples + 1
        self.accuracy_sum = self.accuracy_sum + accuracy
        self.coverage_sum = self.coverage_sum + coverage
        return scores

    def get_accuracy(self):
        self.accuracy_average = self.accuracy_sum / self.total_samples
        return self.accuracy_average
    
    def get_coverage(self):
        self.coverage_average = self.coverage_sum / self.total_samples
        return self.accuracy_average