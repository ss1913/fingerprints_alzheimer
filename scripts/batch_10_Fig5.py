#This code was adapted from https://github.com/gpreti/GSP_StructuralDecouplingIndex/blob/master/Code_NCOMMS/Python/05_metaanalysis_neurosynth_myanalysis.ipynb
#If you use this code please cite the repository and/or associated paper
import os
from neurosynth.base.dataset import Dataset
from neurosynth.analysis import decode
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Set fonttype for SVG
plt.rcParams['svg.fonttype'] = 'none'

# Function to get order
def get_order(data, threshold):
    avg_indices = [np.average(np.array(range(len(item))) + 1, weights=item) for item in data]
    heatmap_order = np.argsort(avg_indices)
    return heatmap_order

# Define data path
data_path = '../functions/data_neurosynth'

# Load dataset
pickled_dataset = os.path.join(data_path, 'dataset.pkl')
dataset = Dataset.load(pickled_dataset)

# Load features
features = pd.read_csv(os.path.join(data_path, 'v3-topics-50.txt'), sep='\t', index_col=0)
topics_to_keep = [ 1, 4,  6, 14, 
                  18, 19, 23, 25, 
                  20, 21, 27, 29,
                  30, 31, 33, 35, 
                  36, 38, 37, 41, 
                  44, 45, 48, 49]
labels = ['face/affective processing', ' verbal semantics', 'cued attention', 'working memory', 
          'autobiographical memory', 'reading', 'inhibition', 'motor', 
          'visual perception', 'numerical cognition', 'reward-based decision making', 'visual attention', 
          'multisensory processing', 'visuospatial','eye movements', 'action',
          'auditory processing', 'pain', 'language', 'declarative memory', 
          'visual semantics', 'emotion', 'cognitive control', 'social cognition']
features = features.iloc[:, topics_to_keep]
features.columns = labels
dataset.add_features(features, append=False)

# Set up decoder
decoder = decode.Decoder(dataset, method='roi')

# Set threshold
threshold = 3.1

# Define paths
nifti_dir = '../outputs/neurosynth'
decoding_paths = [
    ('ICC_sig_GENEVA_ADNI_CU_minus_shen_dil_bin.nii.gz',
     'decoding_results_CU_minus_overlap.txt'),
    ('ICC_sig_GENEVA_ADNI_MCI_plus_shen_dil_bin.nii.gz',
     'decoding_results_MCI_plus_overlap.txt'),
    ('ICC_sig_GENEVA_ADNI_Dementia_plus_shen_dil_bin.nii.gz',
     'decoding_results_Dementia_plus_overlap.txt')
]

# Perform decoding and save results
for input_file, save_file in decoding_paths:
    input_path = os.path.join(nifti_dir, input_file)
    save_path = os.path.join(nifti_dir, save_file)
    data = decoder.decode(input_path, save=save_path)

    # Apply threshold
    data[data < threshold] = 0.01

    # Reorder columns
    heatmap_order = get_order(np.array(data), threshold)
    plot_data = data.reindex(data.index[heatmap_order])

    # Plot heatmap
    sns.set(context="paper", font="sans-serif", font_scale=1)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10), sharey=True)
    sns.heatmap(plot_data, linewidths=1, square=True, cmap='Greys', robust=False,
                ax=ax, vmin=2, vmax=7, mask=plot_data == 0)
    cbar = ax.collections[0].colorbar
    cbar.set_label('z-stat', rotation=270)
    cbar.set_ticks(ticks=[threshold, 7])
    cbar.set_ticklabels(ticklabels=[threshold, 7])
    cbar.outline.set_edgecolor('black')

    # Save figure
    fig.savefig(os.path.join(nifti_dir, save_file.replace('.txt', '.svg')), format='svg')
