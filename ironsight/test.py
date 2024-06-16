from pathlib import Path
from Predict_Ferro_Betas import Predict_Ferro_Betas

# Define the parameters for Meth_Seg
# bam_file_list = ['/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_NB2/IMR21nov_R2.sorted.bam',
#                  '/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_NB2/SHP17nov_R2.sorted.bam',
#                  '/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_20231114/no_haplotag/SHY2sept_R1.sorted.bam',
#                  '/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_20231114/no_haplotag/SKN14nov_R2.sorted.bam']
# ids = ['IMR21nov_R2', 'SHP17nov_R2', 'SHY2sept_R1', 'SKN14nov_R2']
bam_file_list = ['/data/lvisser/modkit/outputs/bam/data_OHMX20230016R_NB2/IMR21nov_R2.sorted.bam']
ids = ['IMR21nov_R2']
outdir = '/data/lvisser/ironsight/examples'
model_type = "NB"

Predict_Ferro_Betas(bam_file_list, ids, outdir, model_type)