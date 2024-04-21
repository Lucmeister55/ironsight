from pathlib import Path
from Predict_Ferro import Predict_Ferro

# Define the parameters for Meth_Seg
bam_file_list = [Path('/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_20231114/no_haplotag/SKN14nov_R2.sorted.bam')]
outdir = '/data/lvisser/ironsight/outputs'

Predict_Ferro(bam_file_list, outdir, verbose = False)