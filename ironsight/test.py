from pathlib import Path
from Predict_Ferro import Predict_Ferro

# Define the parameters for Meth_Seg
bam_file_list = ['/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_20231114/no_haplotag/SKN14nov_R2.sorted.bam']
id = ['SKN14nov_R2']
metadata = "/data/lvisser/ironsight/examples/metadata.csv"
outdir = '/data/lvisser/ironsight/examples'

Predict_Ferro(bam_file_list, id, metadata, outdir)