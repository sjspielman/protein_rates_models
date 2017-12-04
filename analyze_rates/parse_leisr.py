"""
    Parse LEISR inferences.
    Produce a file of rates for each and a single CSV of model fits with logl, AICc, Gamma shape parameter.
"""
from phyphy import *
import os
import sys

fit_output = "model_fits.csv"
fitstring = "dataset,model,type,logl,AICc,GammaShape\n"
mito = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]

    
for type in ["organelle", "enzyme"]:
    
    json_directory = "../rate_inference/" + type + "_data-inference/"

    jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]
    
    for json in jsons:
        print json
        p = Extractor(json_directory + json)
    
        output = json.replace(".json", ".csv")
        p.extract_csv(json_directory + output)
    
        model = p.reveal_fitted_models()[0]
    
        dataset = json.split(".")[0]
        outmodelname = json.split(".")[2]
        if "-G" in json:
            shape = str(p.extract_model_rate_distributions(model)["Gamma distribution shape parameter"])
        else:
            shape = "NA"
 
        gene = json.split(".")[0]          
        if type == "organelle":
            if gene in mito:
                finaltype = "mito"
            else:
                finaltype = "chloro"
        else:
            finaltype = "enzyme"

        thisfit = json.split(".")[0] + "," + outmodelname + "," + finaltype + "," + str(p.extract_model_logl(model)) + "," + str(p.extract_model_aicc(model)) + "," + shape +  "\n"
        fitstring += thisfit

fitstring.strip()
with open(fit_output, "w") as f:
    f.write(fitstring)
    