'''
python script to prepare lists of raw and targets vcf files for input into the variant calling pipeline bcftools merge steps
'''
import glob
import os
import sys
def prepare_vcf_list_files(output_dir):
    # define variant calling dir patterns
    patterns = [
        os.path.join(output_dir, 'variant_calling', '*_raw', 'merge_output.vcf.gz'),
        os.path.join(output_dir, 'variant_calling', '*_targets', 'merge_output.vcf.gz')
    ]
    # make lists to store raw and target filepaths
    vcf_paths_raw = []
    vcf_paths_targets = []
    #get vcf paths for all samples raw and targets vcfs
    for pattern in patterns:
        for vcf_path in glob.glob(pattern):
            if "_raw" in vcf_path:
                vcf_paths_raw.append(vcf_path)
            elif "_targets" in vcf_path:
                vcf_paths_targets.append(vcf_path)
    # write all samples raw and targets paths to two separate txt files
    raw_output_path = os.path.join(output_dir, 'vcf_paths_raw.txt')
    targets_output_path = os.path.join(output_dir, 'vcf_paths_targets.txt')   
    with open(raw_output_path, 'w') as file_raw:
        file_raw.write('\n'.join(vcf_paths_raw))
    with open(targets_output_path, 'w') as file_targets:
        file_targets.write('\n'.join(vcf_paths_targets))
    print(f"VCF lists prepared: {len(vcf_paths_raw)} raw and {len(vcf_paths_targets)} targets")
    print(f"Raw VCF list saved to: {raw_output_path}")
    print(f"Targets VCF list saved to: {targets_output_path}")
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python prepare_vcf_lists.py <output_directory>")
        sys.exit(1)
    output_dir = sys.argv[1]
    prepare_vcf_list_files(output_dir)
