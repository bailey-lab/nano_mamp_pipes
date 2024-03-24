project_resources='/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_project_resources/DR23KE'
species_resources='/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7'
output_data='/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/variant_calling'
sif_file='/nfs/jbailey5/baileyweb/bailey_share/bin/miptools_dev_24-02-13.sif'

# this will be replaced with non-notebook version that automatically runs all necessary code from notebook and outputs the tables
# use vcf_to_tables() in analysis-of-test-data-GATK.ipynb not freebayes version due to differences in AD formatting
# will need to include probe_set specific config options for probe sets with and without gene_name and amino acid change columns in targets file (aggregate amino acids vs nucleotides have different settings in vcf_to_tables())
singularity run \
  -B $project_resources:/opt/project_resources \
  -B $species_resources:/opt/species_resources \
  -B $output_data:/opt/analysis \
  --app jupyter $sif_file

# Put miplithon and mamplicorn commands here or in subsequent script step