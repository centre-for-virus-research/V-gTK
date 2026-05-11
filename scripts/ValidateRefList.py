# check there are three columns, 
#1st is accession, 2nd contains master or reference, 3rd is segment (if segmented)
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Validate reference list format.")
parser.add_argument("-r", "--ref_list", required=True, help="Path to the reference list TSV file.")
parser.add_argument("-s", "--is_segmented", required=True, help="Indicates if the data is segmented ('Y' or 'N').")
parser.add_argument("-o", "--output", required=False, help="Path to the output file.")
params = parser.parse_args()

# read tsv
df = pd.read_csv(params.ref_list, sep="\t", header=None)

if df.shape[1] < 2:
    raise ValueError("Reference list must contain at least two columns: accession and master/reference indicator.")
if params.is_segmented == 'Y' and df.shape[1] < 3:
    raise ValueError("For segmented data, reference list must contain three columns: accession, master/reference indicator, and segment.")

# check master/reference column contains only 'master', 'reference', or 'exclusion_list'
valid_indicators = {'master', 'reference', 'exclusion_list'}
if not df[1].isin(valid_indicators).all():
    invalid = set(df[1].unique()) - valid_indicators
    raise ValueError(f"Indicator column contains invalid values: {invalid}. Allowed: {valid_indicators}")

# check for each unique segment value there's a master (unless it's an exclusion-only segment)
if params.is_segmented == 'Y':
    for segment in df[2].unique():
        segment_df = df[df[2] == segment]
        segment_types = set(segment_df[1].values)
        # Segments that only contain exclusion_list entries don't need a master
        if segment_types == {'exclusion_list'}:
            print(f"Segment '{segment}' is an exclusion-only segment (all entries are exclusion_list). Skipping master check.")
            continue
        if 'master' not in segment_df[1].values:
            raise ValueError(f"No master found for segment '{segment}'. Each non-exclusion segment must have one master entry.")

print("Reference list validation passed.")
if params.output:
    with open(params.output, 'w') as out_file:
        out_file.write("Reference list validation passed.\n")