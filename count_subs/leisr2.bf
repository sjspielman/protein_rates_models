RequireVersion("2.3.5");


LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestralSJS.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/models/DNA.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/models/DNA/HKY85.bf");
LoadFunctionLibrary("libv3/models/DNA/JC69.bf");
LoadFunctionLibrary("libv3/models/protein.bf");
LoadFunctionLibrary("libv3/models/protein/empirical.bf");

// for JSON storage compatibility
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");

/*------------------------------------------------------------------------------*/

utility.ToggleEnvVariable ("NORMALIZE_SEQUENCE_NAMES", 1);

leisr.analysis_description = {
    terms.io.info: "LEISR (Likelihood Estimation of Individual Site Rates) infer relative amino-acid or nucleotide rates from a fixed nucleotide or amino-acid alignment and tree, with possibility for partitions. Relative site-specific substitution rates are
    inferred by first optimizing alignment-wide branch lengths, and then inferring a site-specific uniform tree scaler.",
    terms.io.version: "0.4",
    terms.io.reference: "Spielman, S.J. and Kosakovsky Pond, S.L. (2018). Relative evolutionary rate inference in HyPhy with PeerJ 6:e4339. DOI 10.7717/peerj.4339 ; Pupko, T., Bell, R. E., Mayrose, I., Glaser, F. & Ben-Tal, N. Rate4Site: an algorithmic tool for the identification of functional regions in proteins by surface mapping of evolutionary determinants within their homologues. Bioinformatics 18, S71–S77 (2002).",
    terms.io.authors: "Sergei L Kosakovsky Pond and Stephanie J Spielman",
    terms.io.contact: "{spond,stephanie.spielman}@temple.edu"
};

io.DisplayAnalysisBanner(leisr.analysis_description);
/*******************************************************************************************************************/



/***************************************** LOAD DATASET **********************************************************/
SetDialogPrompt ("Specify a multiple sequence alignment file");
leisr.alignment_info  = alignments.ReadNucleotideDataSet ("leisr.dataset", None);

leisr.name_mapping = leisr.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];
if (None == leisr.name_mapping) {
    leisr.name_mapping = {};
    utility.ForEach (alignments.GetSequenceNames ("leisr.dataset"), "_value_", "`&leisr.name_mapping`[_value_] = _value_");
}


leisr.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (leisr.alignment_info[utility.getGlobalValue("terms.data.partitions")], leisr.name_mapping);
leisr.partition_count = Abs (leisr.partitions_and_trees);
leisr.filter_specification = alignments.DefineFiltersForPartitions (leisr.partitions_and_trees, "leisr.dataset" , "leisr.filter.", leisr.alignment_info);

/*******************************************************************************************************************/


/***************************************** MODEL SELECTION **********************************************************/

leisr.protein_type    = "Protein";
leisr.nucleotide_type = "Nucleotide";
leisr.analysis_type  = io.SelectAnOption ({{leisr.protein_type , "Infer relative rates from a protein (amino-acid) alignment"}, {leisr.nucleotide_type, "Infer relative rates from a nucleotide alignment"}},
                                                    "Select your analysis type:");
                                                    

if (leisr.analysis_type ==  leisr.protein_type) {
    leisr.baseline_model  = io.SelectAnOption (models.protein.empirical_models, "Select a protein model:");
    leisr.generators = models.protein.empirical.plusF_generators;
}
else {

    leisr.baseline_model  = io.SelectAnOption (models.DNA.models, "Select a nucleotide model:");
    leisr.generators = models.DNA.generators;
}

leisr.use_rate_variation = io.SelectAnOption( {{"Gamma", "Use a four-category discrete gamma distribution when optimizing branch lengths."},
                                                    {"GDD", "Use a four-category general discrete distribution when optimizing branch lengths."},
                                                    {"No", "Do not consider rate variation when optimizing branch lengths."}
                                                    },
                                                    "Optimize branch lengths with rate variation?");


function leisr.Baseline.ModelDescription(type){
    def = Call( leisr.generators[leisr.baseline_model], type);
    return def;
}

function leisr.Baseline.ModelDescription.withGamma(type){
    def = leisr.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.Gamma.factory ({terms.rate_variation.bins : 4});
    return def;
}
function leisr.Baseline.ModelDescription.withGDD4(type){
    def = leisr.Baseline.ModelDescription(type);
	def [terms.model.rate_variation] = rate_variation.types.GDD.factory ({terms.rate_variation.bins : 4});
    return def;
}


leisr.baseline_model_name = leisr.baseline_model;
if (leisr.analysis_type ==  leisr.protein_type) {
    leisr.baseline_model_name  = leisr.baseline_model_name + "F";
}
if (leisr.use_rate_variation == "Gamma"){
    leisr.baseline_model_name      = leisr.baseline_model_name + "+4Gamma";
    leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription.withGamma";
}
else {
    if (leisr.use_rate_variation == "GDD"){
        leisr.baseline_model_name      = leisr.baseline_model_name + "+4GDD";
        leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription.withGDD4";
    }
    else {
        leisr.baseline_model_name      = leisr.baseline_model_name;
        leisr.baseline_model_desc      = "leisr.Baseline.ModelDescription";
    }
}


leisr.outfile  = io.PromptUserForString ("Output file name:");

/*******************************************************************************************************************/



/***************************************** INFERENCE **********************************************************/


leisr.trees = utility.Map (leisr.partitions_and_trees, "_value_", "_value_[terms.data.tree]"); // value => value['tree']
leisr.filter_names = utility.Map (leisr.filter_specification, "_value_", "_value_[terms.data.name]"); // value => value['name']
leisr.alignment_wide_MLES = estimators.FitSingleModel_Ext (
                                                          leisr.filter_names,
                                                          leisr.trees,
                                                          leisr.baseline_model_desc,
                                                          None,
                                                          {terms.run_options.retain_lf_object: TRUE});
leisr.ancestors = ancestral.build (leisr.alignment_wide_MLES[terms.likelihood_function], 0, None);
leisr.site_subs = {"COUNTS":{},
                   "TOTAL":{}};
for (i = 0; i < leisr.alignment_info[terms.data.sites]; i+=1){
    leisr.counts = ancestral._substitutionsBySite(leisr.ancestors,i)["COUNTS"];
    (leisr.site_subs["COUNTS"])[i] = leisr.counts;
    (leisr.site_subs["TOTAL"])[i] = (+leisr.counts);

}

io.SpoolJSON (leisr.site_subs, leisr.outfile);



