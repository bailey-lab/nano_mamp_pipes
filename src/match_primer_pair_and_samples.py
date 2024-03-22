import sys
import pandas as pd
import os
# assumed two input files are formatted same as the two files Alec sent George and Alfred on 240229
def main(samples_file, primer_pairs_file,output_dir):
    # read in the primer pair index and sample sheet template
    samples_df = pd.read_csv(samples_file)
    primer_pairs_df = pd.read_csv(primer_pairs_file)
    # merge the dataframes on FW_plate, REV_plate, and well_location
    merged_df = pd.merge(samples_df, primer_pairs_df, on=['FW_plate', 'REV_plate', 'well_location'], how='inner')
    # select results columns
    final_df = merged_df[['sample', 'primer_pair']]
    # define output filepath
    output_file = os.path.join(output_dir, 'matched_samples_primer_pairs.csv')
    # make output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # save results table to csv
    final_df.to_csv(output_file, index=False)
    # success message
    print('Matched samples and primer pairs have been saved to matched_samples_primer_pairs.csv')
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python match_tables.py samples.csv primer_pairs.csv output_directory")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])