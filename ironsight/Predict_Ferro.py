# -*- coding: utf-8 -*-

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
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder

# Local imports
from ironsight import __version__ as version
from ironsight.common import *


# ~~~~~~~~~~~~~~~~~~~~~~~~Helper Functions~~~~~~~~~~~~~~~~~~~~~~~~#

def Predict_Ferro(modbam_list: [str],
    id_list: [str],
    metadata: str,
    outdir: str,
    verbose: bool = False,
    quiet: bool = False,
    progress: bool = False,
    **kwargs,
):
    # Init method
    opt_summary_dict = opt_summary(local_opt=locals())
    log = get_logger(name="ironsight_Predict_Ferro", verbose=verbose, quiet=quiet)
    
    log.warning("Checking options and input files")
    log_dict(opt_summary_dict, log.debug, "Options summary")
    
    # At least one output file is required, otherwise it doesn't make any sense
    log.debug("Checking required output")
    
    if not outdir:
        raise ironsightError("At least 1 output file is requires (-t or -b)")

    # Load your pre-trained model
    model = load('models/model_wgbs_nb1_train_0_1_2_test_3.joblib')

    metadata = pd.read_csv(metadata) 

    # Convert the BAM file to a format that can be used as input to your model
    for i, modbam in enumerate(modbam_list):
        id = id_list[i]
        meth_seg_fm = convert_bam_to_input(modbam, metadata[metadata["id"] == id], model)

    X_test = meth_seg_fm.select_dtypes(exclude=['object'])

    # Make a prediction using the model
    prediction = model.predict(meth_seg_fm)

    # Generate an HTML report
    generate_html_report(prediction)

def convert_bam_to_input(modbam, metadata, model):
    # This function should convert the BAM file to the format expected by your model
    with tempfile.TemporaryDirectory() as tmp_dir:
        print('Temp directory name:', tmp_dir)

    call_bash_script("extract_table.sh", modbam, tmp_dir)

    cpg_cat = process_and_concat_tables(tmp_dir)

    # Assuming `model` is your trained model
    if hasattr(model, 'get_feature_names_out'):
        feature_names = model.get_feature_names_out()

        # Create a DataFrame from feature_names
        segments = pd.DataFrame(feature_names, columns=['segment_id'])

        # Add some additional columns
        segments['chrom'] = segments['segment_id'].str.extract(r'chr(\d+)')
        segments['start'] = segments['segment_id'].str.extract(r':(\d+)-\d+')
        segments['end'] = segments['segment_id'].str.extract(r':\d+-(\d+)')
        segments['gene_symbol'] = None
        segments['length'] = segments['end'].astype(int) - segments['start'].astype(int)
    else:
        print("The model does not support getting feature names.")

    meth_seg = cpg2segment_aggregation_trees(cpg_cat, segments[["chrom", "start", "end", "segment_id", "gene_symbol", "length"]].copy())

    meth_seg_fm = create_fm(meth_seg, metadata)

    return meth_seg_fm

def generate_html_report(prediction):
    # Load the template
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('ironsight/templates/report.html')

    # Render the template with the prediction
    report = template.render(prediction=prediction)

    # Write the report to a file
    with open('ironsight/templates/report.html', 'w') as f:
        f.write(report)

def process_and_concat_tables(input_directory, metadata = None, id = None):
    files = [os.path.join(input_directory, file) for file in os.listdir(input_directory)]
    for file in tqdm(files):
        temp_rawdata = pd.read_csv(file, sep='\t')
        temp_rawdata = temp_rawdata.drop(columns = ["forward_read_position", "fw_soft_clipped_start", "fw_soft_clipped_end", "ref_kmer", "query_kmer", "canonical_base", "modified_primary_base", "inferred", "flag"])
        temp_rawdata = temp_rawdata[temp_rawdata['ref_position'] != -1]
        temp_rawdata = temp_rawdata[temp_rawdata['base_qual'] > 20]
        temp_rawdata = temp_rawdata[temp_rawdata['mod_code'] == "m"]
        if id:
            temp_rawdata['id'] = id
            if metadata:
                temp_rawdata.merge(metadata, on='id')
        # Count the number of times each position occurs
        temp_rawdata['cpg_coverage'] = temp_rawdata.groupby(['ref_position', 'chrom'])['ref_position'].transform('count')

        temp_rawdata = temp_rawdata.groupby(['ref_position', 'chrom', 'cpg_coverage'], as_index=False)['mod_qual'].mean()
        temp_rawdata = temp_rawdata.sort_values(['ref_position', 'chrom', 'cpg_coverage']).reset_index(drop=True)

        cpg_cat = pd.concat([cpg_cat, temp_rawdata])

        return cpg_cat
    
def cpg2segment_aggregation_trees(cpg_df, segment_df):

    dataframes = {'cpg_df': ['source_directory', 'sample_id_adj', 'chrom'],
                'segment_df': ['start', 'end', 'segment_id', 'chrom', 'gene_symbol', 'length']}

    for df_name, columns in dataframes.items():
        for column in columns:
            assert_column(eval(df_name), column, df_name)

    # Filter out chromosomes from cpg_df that aren't in segment_df
    cpg_df = cpg_df[cpg_df['chrom'].isin(segment_df['chrom'].unique())]

    # Create an IntervalTree for each group
    segment_trees = {str(seq_id): IntervalTree(Interval(row.start, row.end, row.segment_id) for row in group.itertuples()) 
                     for seq_id, group in segment_df.groupby('chrom')}
    
    print("initializing meth_seg...")

    # Initialize meth_seg DataFrame
    meth_seg = pd.DataFrame(list(itertools.product(cpg_df["sample_id_adj"].unique(), 
                                                   segment_df["segment_id"].unique())),
                            columns=['sample_id_adj', 'segment_id'])
    
    meth_seg = meth_seg.merge(cpg_df[['sample_id_adj', 'source_directory']].drop_duplicates(), on='sample_id_adj', how='left')

    meth_seg = meth_seg.sort_values(['sample_id_adj']).reset_index(drop=True)
    meth_seg["total_methylation"] = 0.0
    meth_seg["positions"] = [[] for _ in range(len(meth_seg))]
    meth_seg["mod_qual_positions"] = [[] for _ in range(len(meth_seg))]
    meth_seg["num_cpgs"] = 0

    # Initialize a dictionary to store the rows
    meth_seg_dict = {(row.segment_id, row.sample_id_adj, row.source_directory): row for _, row in meth_seg.iterrows()}

    # Group the DataFrame by 'source_directory', 'sample_id_adj', and 'chrom'
    grouped = cpg_df.groupby(['source_directory', 'sample_id_adj', 'chrom'])

    for (source_directory, sample_id_adj, chrom), group in tqdm(grouped, desc="Aggregating"):
        tree = segment_trees[chrom]
        for row in group.itertuples():
            intervals = tree[row.ref_position]
            for interval in intervals:
                # Use the dictionary for lookup and update
                key = (interval.data, sample_id_adj, source_directory)
                meth_seg_row = meth_seg_dict[key]
                meth_seg_row.total_methylation += row.mod_qual
                meth_seg_row.num_cpgs += 1  # Increment the number of CpGs for this segment
                meth_seg_row.positions.append(row.ref_position)
                meth_seg_row.mod_qual_positions.append(row.mod_qual)

    # Convert the dictionary back to a DataFrame
    meth_seg = pd.DataFrame(meth_seg_dict.values())

    # Calculate the average of the values in the 'mod_quals_positions' column
    meth_seg['avg_methylation'] = meth_seg['mod_qual_positions'].apply(lambda x: np.mean(x) if x else 0)

    meth_seg = meth_seg.merge(segment_df[['segment_id', 'gene_symbol', 'length', 'chrom']], on='segment_id', how='left').sort_values("segment_id")

    # Convert the lists to strings
    meth_seg['positions'] = meth_seg['positions'].astype(str)
    meth_seg['mod_qual_positions'] = meth_seg['mod_qual_positions'].astype(str)

    # Now you can safely drop duplicates
    meth_seg = meth_seg.drop_duplicates().reset_index(drop=True)

    return meth_seg

def create_fm(meth_seg, metadata):
    dataframes = {
        'meth_seg': ['segment_id', 'source_directory', 'sample_id_adj', 'avg_methylation'],
        'metadata': ['source_directory', 'sample_id_adj', 'Group']
    }

    for df_name, columns in dataframes.items():
        for column in columns:
            assert_column(eval(df_name), column, df_name)

    meth_seg_pivot = meth_seg.pivot_table(index=["sample_id_adj", "source_directory"], columns="segment_id", values="avg_methylation").reset_index()
    # Fill missing values with zeros
    meth_seg_pivot = meth_seg_pivot.fillna(0)
    meth_seg_fm = meth_seg_pivot.merge(metadata[["sample_id_adj", "Group", "tumor_type"]], on = "sample_id_adj", how = "left").drop_duplicates().reset_index(drop = True)
    # Assuming 'Group' is the column you want to move to the front
    group = meth_seg_fm.pop('Group')
    meth_seg_fm.insert(0, 'Group', group)
    meth_seg_fm.sort_values("source_directory", inplace = True)

    # Remove all zero columns
    meth_seg_fm = meth_seg_fm.loc[:, (meth_seg_fm != 0).any(axis=0)]

    return meth_seg_fm