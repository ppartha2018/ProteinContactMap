from __future__ import division
from Bio.PDB import *

class ContactMapMetrics:

    total_samples = 0
    precision_sum = 0
    coverage_sum = 0
    precision_average = 0
    coverage_average = 0

    def cmap_key_equals(self,key1,key2):
    #Key1, Key2 are tuples where each element is of type Bio.PDB.Atom.Atom
        return ((key1[0].__eq__(key2[0])) and (key1[1].__eq__(key2[1])))
    
    '''
    Contact map evaluation measures:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5820169/
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3226919/
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1404-z
    ''' 
    def compare_contact_map(self, computed_map = {}, predicted_map = {}, no_of_atoms = 1, no_of_native_contacts = 1):
        #example
        #computed_map = {('ca1','ca2'):0,('ca1','ca1'):1,('ca1','ca3'):1}
        #predicted_map = {('ca1','ca2'):1,('ca1','ca1'):1,('ca1','ca3'):1}
        
        #print("Length of computed map: %d, Length of predicted map: %d" % (len(computed_map), len(predicted_map)))
        #print computed_map

        precision = 0
        recall = 0
        f1 = 0
        coverage = 0
        
        predicted_map_1s = 0
        predicted_map_0s = 0
        true_positives = 0
        false_positives = 0
        true_negatives = 0
        false_negatives = 0

        for key,value in predicted_map.iteritems():
            #this check: in case of native models, sometimes residues are missing, so don't proceed if such a key does not exist in native computed contact map
            #if key in computed_map:
            if(value == 1):
                predicted_map_1s = predicted_map_1s + 1
                if(computed_map[key] == value):
                    true_positives = true_positives + 1
                else:
                    false_positives = false_positives + 1
            else: # predicated value is 0 - or no prediction
                predicted_map_0s = predicted_map_0s + 1
                if(computed_map[key] == value):
                    true_negatives = true_negatives + 1
                else:
                    false_negatives = false_negatives + 1

        #print computed_map_0s, computed_map_1s, predicted_map_0s, predicted_map_1s, difference
        #precision Acc = TP/(TP + FP), 
        #where TP and FP are the numbers of correctly and incorrectly predicted contacts, respectively
        precision = true_positives / (true_positives + false_positives)
        recall = true_positives / (true_positives + false_negatives)
        print ("TP=%d, FP=%d, TN=%d, FN=%d" % (true_positives,false_positives,true_negatives,false_negatives))
        print ("Evaluation: precision: %f, recall=%f" % (precision,recall))
        try:
            f1 = 2 * precision * recall / (precision + recall)
        except ZeroDivisionError:
            f1 = 0 
        #coverage = TP/ Nc, Nc - no of contacts in native model
        coverage = true_positives / no_of_native_contacts # equal to (TP) / n
        #score2 = (predicted_map_0s + predicted_map_1s) / no_of_atoms
        
        scores = []
        scores.append(precision)
        scores.append(coverage)
        scores.append(recall)
        scores.append(f1)
        #scores.append(score2)
        
        #print ("Evaluation: Score1: %f, Score2: %f, Score3: %f" % (score1,score2,score3))
        print "TP,FN,FP,TN" , true_positives, false_negatives, false_positives, true_negatives
        print ("Evaluation: precision: %f" % (precision))
        self.total_samples = self.total_samples + 1
        self.precision_sum = self.precision_sum + precision
        self.coverage_sum = self.coverage_sum + coverage
        return scores

    def get_precision(self):
        self.precision_average = self.precision_sum / self.total_samples
        return self.precision_average
    
    def get_coverage(self):
        self.coverage_average = self.coverage_sum / self.total_samples
        return self.coverage_average