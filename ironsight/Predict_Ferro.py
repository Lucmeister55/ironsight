# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
import sys
import pandas as pd
from joblib import dump, load
from jinja2 import Environment, FileSystemLoader

# Local imports
from ironsight import __version__ as version
from ironsight.common import *

# ~~~~~~~~~~~~~~~~~~~~~~~~Helper Functions~~~~~~~~~~~~~~~~~~~~~~~~#

def Predict_Ferro(modbam_list: [str],
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
    model = load('ironsight/models/model_wgbs_nb1_train_0_1_2_test_3.joblib') 

    # Convert the BAM file to a format that can be used as input to your model
    input_data = convert_bam_to_input(modbam_list)

    # Make a prediction using the model
    prediction = model.predict(input_data)

    # Generate an HTML report
    generate_html_report(prediction)

def convert_bam_to_input(modbam_list):
    # This function should convert the BAM file to the format expected by your model
    # This is just a placeholder and should be replaced with your actual implementation
    return pd.read_csv(modbam_list)

def generate_html_report(prediction):
    # Load the template
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('ironsight/templates/report.html')

    # Render the template with the prediction
    report = template.render(prediction=prediction)

    # Write the report to a file
    with open('ironsight/templates/report.html', 'w') as f:
        f.write(report)