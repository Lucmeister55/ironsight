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
import pkg_resources
from tqdm import tqdm

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
        for feature in tqdm(feature_names, desc="Writing features to BED file"):
            chrom, positions = feature.split(':')
            chromStart, chromEnd = positions.split('-')
            f.write(f'{chrom}\t{chromStart}\t{chromEnd}\n')

def calculate_beta_values(bam_file, bed_file, output_file, prob_threshold=0.67, mapq_threshold=10):
    try:
        # Load the BED file
        bed = pybedtools.BedTool(bed_file)
        print("Loaded BED file.")

        # Open the output file for writing
        with open(output_file, 'w') as f:
            print(f"Opened output file: {output_file}")

            # Write the header
            f.write('chrom\tchromStart\tchromEnd\tbeta\n')
            print("Wrote header to output file.")

            # Load the BAM file
            bam = pysam.AlignmentFile(bam_file)
            print("Loaded BAM file.")

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
                
                # print(total_cpgs)

                cpg_info = sorted(cpg_info, key=lambda x: x['probability'])

                # Calculate the beta value
                if total_cpgs > 0:
                    beta = high_prob_cpgs / total_cpgs
                    beta = round(beta, 2)
                else:
                    beta = 'NA'

                # Write the beta value to the output file
                f.write(f'{chrom}\t{chromStart}\t{chromEnd}\t{beta}\n')

            print("Finished processing regions.")
    except Exception as e:
        print(f"An error occurred: {e}")

def load_model(model_type, pkg_dir):
    # Load the model
    if model_type == 'NB':
        print("Loading model...")
        model_folder = os.path.join(pkg_dir, 'models/NB/')
        # Get a list of all joblib files in the folder
        joblib_files = glob.glob(model_folder + '*.joblib')
        # Check if there is exactly one joblib file
        if len(joblib_files) != 1:
            raise ValueError('Expected exactly one joblib file in the folder, but found ' + str(len(joblib_files)))
        # Load the joblib file
        model = joblib.load(joblib_files[0])
        print("Model loaded.")
        
        print("Loading explainer...")
        pkl_files = glob.glob(model_folder + '*.pkl')

        # Check if there is exactly one pickle file
        if len(pkl_files) != 1:
            raise ValueError('Expected exactly one pickle file in the folder, but found ' + str(len(pkl_files)))

        # Load the pickle file
        with open(pkl_files[0], 'rb') as f:
            explainer = pickle.load(f)
        print("Explainer loaded.")
    return model, explainer

def calculate_all_beta_values(bam_file_list, ids, temp_dir, bed_file):
    print("Calculating beta values...")
    print(f"bam_file_list: {bam_file_list}")
    print(f"ids: {ids}")
    print(f"temp_dir: {temp_dir}")
    print(f"bed_file: {bed_file}")
    print()
    try:
        # Calculate beta values for each BAM file
        for bam_file, id in zip(bam_file_list, ids):
            print(f"Calculating beta values for {bam_file}...")
            output_file = os.path.join(temp_dir, f'{id}.betas.tsv')
            calculate_beta_values(bam_file, bed_file, output_file)
    except Exception as e:
        print(f"An error occurred: {e}")  # Print any exceptions

def merge_beta_values_with_bed(bed_df, ids, temp_dir):
    # For each id, merge the beta column from the corresponding betas file with the bed file
    for id in tqdm(ids, desc="Merging beta values"):
        betas_df = pd.read_csv(os.path.join(temp_dir, f'{id}.betas.tsv'), sep='\t')
        # Rename the columns of betas_df before the merge operation
        betas_df = betas_df.rename(columns={'beta': f'{id}'})

        # Merge the DataFrames
        bed_df = pd.merge(bed_df, betas_df[['chrom', 'chromStart', 'chromEnd', f'{id}']], on=['chrom', 'chromStart', 'chromEnd'], how='left')

    # Save the merged DataFrame
    bed_df.to_csv(os.path.join(temp_dir, 'betas.tsv'), sep='\t', index=False)
    
    return bed_df

def create_feature_matrix(bed_df, ids):
    print("Creating feature matrix...")
    bed_df["segment_id"] = bed_df.apply(lambda x: f"{x['chrom']}:{x['chromStart']}-{x['chromEnd']}", axis=1)

    feature_matrix = bed_df.drop(columns=["chrom", "chromStart", "chromEnd"])

    print("Pivoting rows to columns...")
    # Pivot rows to columns
    feature_matrix = feature_matrix.pivot_table(columns='segment_id', values=[id for id in ids], dropna=False).reset_index()
    feature_matrix.rename(columns={'index':'id'}, inplace = True)
    
    print("Replacing NA values with 0...")
    # Replace NA values with 0
    feature_matrix = feature_matrix.fillna(0)

    print("Feature matrix created.")
    return feature_matrix

def run_prediction_and_shap_analysis(feature_matrix, model, explainer, results_dir, ids):
    print("Separating features from target variable...")
    # Separate the features from the target variable
    X = feature_matrix.select_dtypes(include=[np.number])

    feature_names = model.feature_names_in_
    X = X[feature_names]

    print("Running prediction...")
    # Run the prediction
    predictions = model.predict(X)

    # Get the prediction probabilities
    probabilities = model.predict_proba(X)

    # Create labels based on the predictions
    labels = ['Resistant' if prediction == 1 else 'Sensitive' for prediction in predictions]

    shap_values = explainer.shap_values(X)

    # Convert SHAP values to a list of lists
    shap_values_list = shap_values.tolist()

    print("Saving predictions, labels, probabilities, and SHAP values to a file...")
    # Save the predictions, labels, probabilities, and SHAP values to a single file

    results = pd.DataFrame({
        'id': ids,
        'prediction': predictions,
        'label': labels,
        'probability_S': probabilities[0].round(2),  # Assuming you want the max probability
        'probability_R': probabilities[1].round(2),  # Assuming you want the max probability
        'shap_values': shap_values_list,  # Add SHAP values as a new column
    })

    results_csv_path = os.path.join(results_dir, 'results.csv')

    results.to_csv(results_csv_path, index=False)

    # Return the predictions, labels, probabilities, and SHAP values
    return results_csv_path, explainer

def predict(bam_file_list, ids, model_type, progress_lines_path, progress_percentage_path, results_dir, temp_dir, pkg_dir):

    steps = [
        'Checking if the BAM files exist and are accessible...',
        'Loading the prediction model...',
        'Writing features to the BED file...',
        'Calculating beta values for the features...',
        'Merging beta values...',
        'Creating the feature matrix...',
        'Running prediction and SHAP analysis...'
    ]

    total_steps = len(steps)
    current_step = 0

    def write_progress():
        nonlocal current_step
        with open(progress_lines_path, 'a') as progress_file:
            progress_file.write(steps[current_step] + '\n')
        with open(progress_percentage_path, 'w') as progress_file:
            progress_file.write(f'{((current_step+1)/total_steps)*100}')
        current_step += 1

    write_progress()

    # Check if the bam files exist and are accessible
    for bam_file in bam_file_list:
        if not os.path.isfile(bam_file):
            raise FileNotFoundError(f"The bam file {bam_file} does not exist or is not accessible.")

    # Remove all files in the temp directory
    files = glob.glob(os.path.join(temp_dir, '*'))
    for f in files:
        os.remove(f)

    # Remove all files in the temp directory
    files = glob.glob(os.path.join(results_dir, '*'))
    for f in files:
        os.remove(f)

    write_progress()
    model, explainer = load_model(model_type, pkg_dir)

    bed_file = os.path.join(temp_dir, 'features.bed')

    write_progress()
    model_features_to_bed(model, bed_file)

    write_progress()
    calculate_all_beta_values(bam_file_list, ids, temp_dir, bed_file)

    # Read the bed file into a DataFrame
    bed_df = pd.read_csv(bed_file, sep='\t', header = None)
    bed_df.columns = ['chrom', 'chromStart', 'chromEnd']

    write_progress()
    bed_df = merge_beta_values_with_bed(bed_df, ids, temp_dir)

    write_progress()
    feature_matrix = create_feature_matrix(bed_df, ids)
    
    write_progress()
    results_csv_path, explainer = run_prediction_and_shap_analysis(feature_matrix, model, explainer, results_dir, ids)

    return results_csv_path, explainer