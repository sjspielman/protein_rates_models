from phyphy import *
import os
import sys
    

for type in ["enzyme"]:

    json_directory = type + "_data-inference/"

    jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]
    
    for jsonf in jsons:
        
        print jsonf
        p = Extractor(json_directory + jsonf)
    
        output = jsonf.replace(".json", ".csv")
        p.extract_csv(json_directory + output)
    

