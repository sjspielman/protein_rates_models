"""
    Parse LEISR inferences.
    Produce a file of rates for each and a single CSV of model fits with logl, AICc, Gamma shape parameter.
"""
from phyphy import *
import os
import sys

output_directory = "summarized_data/"
fit_output = output_directory + "model_fits.csv"
fitstring = "dataset,sequences,sites,model,type,logl,AICc,GammaShape\n"

mito = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
    

for type in ["enzyme", "virus", "organelle"]:

    json_directory = "../rate_inference/" + type + "_data-inference/"

    jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]
    
    for jsonf in jsons:
        
        print jsonf
        p = Extractor(json_directory + jsonf)
    
        output = jsonf.replace(".json", ".csv")
        p.extract_csv(json_directory + output)
    
        sites = str(p.extract_number_sites())
        seqs = str(p.extract_number_sequences())
           
        jsonmodel = p.reveal_fitted_models()[0]

        dataset = jsonf.split(".")[0]
        outmodelname = jsonf.split(".")[2]
        if "-G" in jsonf:
            shape = str(p.extract_model_rate_distributions(jsonmodel)["Gamma distribution shape parameter"])
        else:
            shape = "NA"
 
        gene = jsonf.split(".")[0]          
        if type == "organelle":
            if gene in mito:
                finaltype = "mito"
            else:
                finaltype = "chloro"
        else:
            finaltype = type

        thisfit = dataset + "," + seqs + "," + sites + "," + outmodelname + "," + finaltype + "," + str(p.extract_model_logl(jsonmodel)) + "," + str(p.extract_model_aicc(jsonmodel)) + "," + shape +  "\n"
        fitstring += thisfit

fitstring.strip()
with open(fit_output, "w") as f:
    f.write(fitstring)
    
    
    
########################################################################################

######### The custom matrix data ############

output_directory = "summarized_data/"
fit_output = output_directory + "model_fits_custom.csv"
fitstring = "dataset,sequences,sites,model,logl,AICc\n"

json_directory = "../rate_inference/custom_matrix-inference/"
jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]

for jsonf in jsons:
    
    print jsonf
    p = Extractor(json_directory + jsonf)

    output = jsonf.replace(".json", ".csv")
    p.extract_csv(json_directory + output)

    sites = str(p.extract_number_sites())
    seqs = str(p.extract_number_sequences())
       
    jsonmodel = p.reveal_fitted_models()[0]

    dataset = jsonf.split(".")[0]
    outmodelname = jsonf.split(".")[2]

    gene = jsonf.split(".")[0]          

    thisfit = dataset + "," + seqs + "," + sites + "," + outmodelname  + "," + str(p.extract_model_logl(jsonmodel)) + "," + str(p.extract_model_aicc(jsonmodel)) +  "\n"
    fitstring += thisfit

fitstring.strip()
with open(fit_output, "w") as f:
    f.write(fitstring)

    