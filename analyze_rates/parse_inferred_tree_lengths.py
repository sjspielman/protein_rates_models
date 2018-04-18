from phyphy import *
import os
import sys

output_directory = "summarized_data/"
lines = "dataset,model,rv,type,treelength\n"

mito = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]

for type in ["gpcr", "enzyme", "organelle"]:

    json_directory = "../rate_inference/" + type + "_data-inference/"
    jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]
    
    for jsonf in jsons:
        
        print jsonf
        p = Extractor(json_directory + jsonf)
        
        jsonmodel = p.reveal_fitted_models()[0]
        dataset = jsonf.split(".")[0]
        modelname = jsonf.split(".")[2]
        model = modelname.split("-")[0]
        rv = modelname.split("-")[1]
        
        bld = p.extract_branch_attribute(jsonmodel)
        tl = 0
        for bl in bld:
            tl += float(bld[bl])

        gene = jsonf.split(".")[0]          
        if type == "organelle":
            if gene in mito:
                finaltype = "mito"
            else:
                finaltype = "chloro"
        else:
            finaltype = type

        thisline = dataset + "," + model + "," + rv + "," + finaltype + "," + str(tl) +  "\n"
        lines += thisline

lines.strip()
with open(output_directory + "inferred_tree_lengths.csv", "w") as f:
    f.write(lines)
    
    
