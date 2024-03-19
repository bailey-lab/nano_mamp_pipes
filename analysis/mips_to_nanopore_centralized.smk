# Need to reformat the rule_all() in each sub smk so that it is conditionally defined if yaml has definition of SANDALONE=TRUE. 
# This is needed if using import smk1 smk2 smk3 because can not have identically named rule_all defined at once. 
# Only define rule all in central smk if STANDALONE=FALSE. Run the individual SMK.
# Need to consider if we will run from central SMK even when STANDALONE=TRUE. This way users dont need to both modify the yaml STANDALONE=FALSE definition and also run a different smk than usual. Maybe this is handled by singularity.
