"""
    Output CSV
        model,site,subType,Count
"""
import os
import json
import numpy as np


aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
output_directory = "../analyze_rates/summarized_data/"


mito = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
    

for mytype in ["gpcr", "enzyme", "organelle"]:

    outfile = output_directory + mytype + "_substitution_counts.csv"
    outstring = "dataset,type,model,site,subtype,count\n"

    raw_directory = mytype + "_data-inference/"
    files = [x for x in os.listdir(raw_directory) if x.endswith("txt")]
    for file in files:
        print file
        dataset = file.split(".")[0]
        model = file.split(".")[2].split("-")[0]
        if mytype == "organelle":
            if dataset in mito:
                finaltype = "mito"
            else:
                finaltype = "chloro"
        else:
            finaltype = mytype

        with open(raw_directory + file, "r") as f:
            raw = json.load(f)
        
        ######## Parse total #########
        total_dict = raw["TOTAL"]
        
        for i in range(len(total_dict)):
            outstring += ",".join([dataset, finaltype, model, str(i+1), "total", str(total_dict[str(i)])]) + "\n"
        
        ####### Parse specific substitutions,saving only where total > 1 #######
        count_dict = raw["COUNTS"]
        
        for i in range(len(count_dict)):
            site_dict = count_dict[str(i)]
            ## we need to do full matrix, not just lower triangle
            for x in range(20):
                for y in range(x+1,20):
                    pair = aminos[x] + aminos[y] 
                    paircount = site_dict[x][y] + site_dict[y][x]
                    if paircount == 0:
                        continue
                    else:
                        outstring += ",".join([dataset, finaltype, model, str(i+1), pair, str(paircount)]) + "\n"
                    
    with open(outfile, "w") as f:
        f.write(outstring.strip())