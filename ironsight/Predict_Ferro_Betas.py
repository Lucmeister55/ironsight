# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
import sys
import pandas as pd
import numpy as np
import shap
import multiprocessing as mp
import seaborn as sns
import matplotlib.pyplot as plt
from joblib import dump, load
from jinja2 import Environment, FileSystemLoader
import tempfile
import os
import re
import pandas as pd
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
import itertools
import sklearn
import joblib
import pysam
import pybedtools
import glob
import shap
import pickle
import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

# Local imports
from ironsight import __version__ as version
from ironsight.common import *

# ~~~~~~~~~~~~~~~~~~~~~~~~Helper Functions~~~~~~~~~~~~~~~~~~~~~~~~#

def model_features_to_bed(model, bed_file):

    # Extract feature names
    feature_names = model.feature_names_in_

    # Open the BED file for writing
    with open(bed_file, 'w') as f:
        # Write each feature name to the BED file
        for feature in feature_names:
            chrom, positions = feature.split(':')
            chromStart, chromEnd = positions.split('-')
            f.write(f'{chrom}\t{chromStart}\t{chromEnd}\n')

def calculate_beta_values(bam_file, bed_file, output_file, prob_threshold=0.67, mapq_threshold=10):
    # Load the BED file
    bed = pybedtools.BedTool(bed_file)

    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Write the header
        f.write('chrom\tchromStart\tchromEnd\tbeta\n')

        # Load the BAM file
        bam = pysam.AlignmentFile(bam_file)

        # Intersect the BAM file with the BED file
        for region in bed:
            chrom, chromStart, chromEnd = region.fields[:3]
            chromStart = int(chromStart)
            chromEnd = int(chromEnd)
            total_cpgs = 0
            high_prob_cpgs = 0

            # Initialize an empty list to store the read IDs
            reads = []

            cpg_info = []

            EXCLUDE_FLAGS = 1796  # Replace this with the sum of the flags you want to exclude

            # Process each read in the region
            for read in bam.fetch(chrom, chromStart, chromEnd):
                # Add the read ID to the list
                reads.append(read)

                # Skip unmapped reads and reads with low mapping quality
                if read.mapping_quality <= mapq_threshold or (read.flag & EXCLUDE_FLAGS):
                    continue

                # Calculate the number of CpGs in the read and the number of CpGs with high probability
                try:
                    # Get the reference positions of the read
                    ref_positions = read.get_reference_positions(full_length = True)

                    ref_positions = [pos + 1 if pos is not None else None for pos in ref_positions]

                    # Get the strand of the read
                    strand = 1 if read.is_reverse else 0

                    if read.modified_bases:
                        try:
                            modified_bases = [(ref_positions[pos], round(prob / 255, 2)) for pos, prob in read.modified_bases[('C', strand, 'm')] if pos < len(ref_positions)]
                        except IndexError:
                            print(f"Warning: MM tag refers to bases beyond sequence length for read {read.query_name}")
                            continue

                        # ref_seq = read.get_reference_sequence()
                        quality_scores = read.query_qualities

                        for i, m in enumerate(modified_bases):
                            if m[0] is not None and chromStart <= m[0] <= chromEnd:
                                total_cpgs += 1
                                index = ref_positions.index(m[0])
                                # Save the information about the position
                                cpg_info.append({
                                    'position': m[0],
                                    'probability': m[1],
                                    'read_id': read.query_name,
                                    'quality_score': quality_scores[index],
                                })
                                if m[1] > prob_threshold:
                                    high_prob_cpgs += 1
                except KeyError:
                    continue
            
            print(total_cpgs)

            cpg_info = sorted(cpg_info, key=lambda x: x['probability'])

            # Calculate the beta value
            if total_cpgs > 0:
                beta = high_prob_cpgs / total_cpgs
                beta = round(beta, 2)
            else:
                beta = 'NA'

            # Write the beta value to the output file
            f.write(f'{chrom}\t{chromStart}\t{chromEnd}\t{beta}\n')

def Predict_Ferro_Betas(bam_file_list, ids, outdir, model_type):

    # Get the directory of the current file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Create the output directory if it does not exist
    os.makedirs(outdir, exist_ok=True)

    # Load the model
    if model_type == 'NB':
        model_folder = os.path.join(current_dir, 'models/NB/')
        # Get a list of all joblib files in the folder
        joblib_files = glob.glob(model_folder + '*.joblib')
        # Check if there is exactly one joblib file
        if len(joblib_files) != 1:
            raise ValueError('Expected exactly one joblib file in the folder, but found ' + str(len(joblib_files)))
        # Load the joblib file
        model = joblib.load(joblib_files[0])
        pkl_files = glob.glob(model_folder + '*.pkl')

        # Check if there is exactly one pickle file
        if len(pkl_files) != 1:
            raise ValueError('Expected exactly one pickle file in the folder, but found ' + str(len(pkl_files)))

        # Load the pickle file
        with open(pkl_files[0], 'rb') as f:
            explainer = pickle.load(f)

    bed_file = os.path.join(current_dir, 'feature_matrices/features.bed')
    model_features_to_bed(model, bed_file)

    # Calculate beta values for each BAM file
    for bam_file, id in zip(bam_file_list, ids):
        output_file = os.path.join(current_dir, f'feature_matrices/{id}.betas.tsv')
        calculate_beta_values(bam_file, bed_file, output_file)

    # Read the bed file into a DataFrame
    bed_df = pd.read_csv(bed_file, sep='\t', header = None)
    bed_df.columns = ['chrom', 'chromStart', 'chromEnd']

    # For each id, merge the beta column from the corresponding betas file with the bed file
    for id in ids:
        betas_df = pd.read_csv(os.path.join(current_dir, f'feature_matrices/{id}.betas.tsv'), sep='\t')
        # Rename the columns of betas_df before the merge operation
        betas_df = betas_df.rename(columns={'beta': f'{id}'})

        # Merge the DataFrames
        bed_df = pd.merge(bed_df, betas_df[['chrom', 'chromStart', 'chromEnd', f'{id}']], on=['chrom', 'chromStart', 'chromEnd'], how='left')

    # Save the merged DataFrame
    bed_df.to_csv(os.path.join(current_dir, 'feature_matrices/betas.tsv'), sep='\t', index=False)

    bed_df["segment_id"] = bed_df.apply(lambda x: f"{x['chrom']}:{x['chromStart']}-{x['chromEnd']}", axis=1)

    feature_matrix = bed_df.drop(columns=["chrom", "chromStart", "chromEnd"])

    # Pivot rows to columns
    feature_matrix = feature_matrix.pivot_table(columns='segment_id', values=[id for id in ids], dropna=False).reset_index()
    feature_matrix.rename(columns={'index':'id'}, inplace = True)
    # Replace NA values with 0
    feature_matrix = feature_matrix.fillna(0)

    # Separate the features from the target variable
    X = feature_matrix.select_dtypes(include=[np.number])

    feature_names = model.feature_names_in_
    X = X[feature_names]

    # Run the prediction
    predictions = model.predict(X)

    # Compute the SHAP values
    shap_values = explainer.shap_values(X)

    # Create a SHAP summary plot
    shap.summary_plot(shap_values, X, feature_names = X.columns, show=False, plot_size = [7, 4])
    plt.savefig(os.path.join(outdir, 'summary_plot.png'))
    plt.clf()  # Clear the current figure

    # Create a SHAP force plot for the first prediction
    shap.plots.force(explainer.expected_value, shap_values[0])
    plt.savefig(os.path.join(outdir, 'force_plot_0.png'))
    plt.clf()  # Clear the current figure

    print(predictions)