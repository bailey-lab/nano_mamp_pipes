'''
python script to prepare lists of raw and targets gvcf files for input into the variant calling pipeline bcftools merge steps
'''
import glob
import os
import sys
def prepare_gvcf_list_files(output_dir):
    # define variant calling dir patterns
    patterns = [
        os.path.join(output_dir, 'variant_calling', '*_raw', 'merge_output.gvcf.gz'),
        os.path.join(output_dir, 'variant_calling', '*_targets', 'merge_output.gvcf.gz')
    ]
    # make lists to store raw and target filepaths
    gvcf_paths_raw = []
    gvcf_paths_targets = []
    #get gvcf paths for all samples raw and targets gvcfs
    for pattern in patterns:
        for gvcf_path in glob.glob(pattern):
            if "_raw" in gvcf_path:
                gvcf_paths_raw.append(gvcf_path)
            elif "_targets" in gvcf_path:
                gvcf_paths_targets.append(gvcf_path)
    # write all samples raw and targets paths to two separate txt files
    raw_output_path = os.path.join(output_dir, 'gvcf_paths_raw.txt')
    targets_output_path = os.path.join(output_dir, 'gvcf_paths_targets.txt')   
    with open(raw_output_path, 'w') as file_raw:
        file_raw.write('\n'.join(gvcf_paths_raw))
    with open(targets_output_path, 'w') as file_targets:
        file_targets.write('\n'.join(gvcf_paths_targets))
    print(f"gvcf lists prepared: {len(gvcf_paths_raw)} raw and {len(gvcf_paths_targets)} targets")
    print(f"Raw gvcf list saved to: {raw_output_path}")
    print(f"Targets gvcf list saved to: {targets_output_path}")
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python prepare_gvcf_lists.py <output_directory>")
        sys.exit(1)
    output_dir = sys.argv[1]
    prepare_gvcf_list_files(output_dir)

# for testing - define directory to input into the function (the root output directory from the pipeline)
# output_dir = "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output"
# run it
# prepare_gvcf_list_files(output_dir)
