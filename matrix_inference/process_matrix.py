import re
import os    
import json
import numpy as np
import pprint
import sys

ZERO=1e-8
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

#####
class JSON_terms():
    def __init__(self):
        
        self.global_          = "global"
        self.logl             = "LogL"
        self.rev_phase_prefix = "REV-Phase-"
        self.final_phase      = "REV-Final"
        self.MLE              = "MLE"
        self.tree             = "Trees"
        self.aa_rate_regex    = u"Substitution rate from amino-acid (\w) to amino-acid (\w)" ### JSONs in unicode
        self.efv_regex        = u"Equilibrium frequency for (\w)"
        self.EFV              = "EFV"
        
    def aa_rate(self, source, target):
        return "Substitution rate from amino-acid " + source + " to amino-acid " + target
        
    def freq(self, AA):
        return "Equilibrium frequency for " + str(AA)
        
        
        
    
def rev_json_to_array( global_d, terms ):
    """
        Parse REV matrix from its JSON structure into a numpy array.
    """
    # Empty matrix, also ensures that diagonal will be 0
    rev = np.zeros([20,20])
  
    for key in global_d:
        if key.startswith("Substitution"):
            find_to_from = re.search(terms.aa_rate_regex, key)
            source = None
            target = None
            if find_to_from:
                source = AMINO_ACIDS.index( find_to_from.group(1))
                target = AMINO_ACIDS.index( find_to_from.group(2))
            assert(source is not None and target is not None), "\nCould not parse a REV rate from JSON dictionary."
        
            try:
                rate = float( global_d[key][terms.MLE] )
            except:
                raise KeyError, "\nCould not parse a REV rate from JSON dictionary."
        
            rev[source][target] = rate
            rev[target][source] = rate   
 
    assert(np.array_equal(rev.diagonal(),np.zeros(20))), "\nDiagonal of a parsed matrix is not 0"
    
    return rev 




def efv_json_to_array( global_d, terms ):
    """
        Parse EFV vector from its JSON structure into a numpy array.
    """

    # Empty matrix, also ensures that diagonal will be 0
    efv = np.zeros(20)
  
    for key in global_d:
        if key.startswith("Equilibrium"):
            find_efv = re.search(terms.efv_regex, key)
            AA = None
            if find_efv:
                AA = AMINO_ACIDS.index( find_efv.group(1))
            assert(AA is not None), "\nCould not parse an EFV from JSON."
            
            try:
                f = float( global_d[key][terms.MLE] )
            except:
                raise KeyError, "\nCould not parse an EFV rate from JSON."
            efv[AA] = f
            
    return efv



def matrix_to_hyphy(matrix):

    rd = {}

    for s in range(len(AMINO_ACIDS)):
        d = {}
        for t in range(s+1, len(AMINO_ACIDS)):
            d[AMINO_ACIDS[t]] = matrix[s][t]
        rd[AMINO_ACIDS[s]] = d
    pprint.pprint(rd)
    



def main():

    terms = JSON_terms()
     
    json_file = sys.argv[1]  #"rpoC1.list.json"        
    with open(json_file, "r") as f:
        parsed_json = json.load(f)
    final_phase = parsed_json[terms.final_phase][terms.global_]
    
    matrix = rev_json_to_array(final_phase, terms)
    #frequencies = efv_json_to_array(final_phase, terms)
    
    ### Add this printed business to hyphy, no need for frequencies because they will be the same by definition.. ###
    matrix_to_hyphy(matrix)
    
#     for headkey in headkeys:
#             if str(headkey).startswith("REV-"):   # or str(headkey) == "WAG-Phase-0":
#             
#                 phase = parsed_json[ headkey ]
# 
#                 ## Q matrix ##
#                 q = rev_json_to_array ( phase[terms.global_], terms )
#                 Q_matrices[ str(headkey) ] = q
#                 Q_rmse = calculate_rmse(q, LG_MATRIX)
#                 logL = float( phase[terms.logl] )
# 
# 
# "REV-Final":{
#    "global":{
#      "Equilibrium frequency for A":{
#        "ID":"vofaSMEV.model.equilibrium_frequency_of.A",
#        "MLE":0.04379282652822862
#       },
# 
# 
# 

main()
