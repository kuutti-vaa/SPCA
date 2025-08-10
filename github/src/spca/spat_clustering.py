import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from esda.moran import Moran
from pysal.lib import weights
from scipy.spatial.distance import cdist
from .data_proc import img_conversion

def pysal_w_conv(weight_array, threshold = False):
    """
    converts a weight matrix np array to pysal weight object. weight array is of shape (num points, num_points)
    Args:
    weight_array (numpy array): numpt array of weights
    threshold: weight threshold below which weight is negligible
    """
    n_points = weight_array.shape[0]
    # Create neighbors and weights dicts
    neighbors_dict = {}
    weights_dict = {}

    for i in range(n_points):
        # Get nonzero neighbors (could filter small weights if needed)
        if threshold: neighbors = np.where(weight_array[i] > threshold)[0]
        else: neighbors = np.where(weight_array[i] > 0)[0]
        weights_i = weight_array[i, neighbors]

        neighbors_dict[i] = list(neighbors)
        weights_dict[i] = list(weights_i)

    # Now create a W object
    w_gaussian = weights.W(neighbors_dict, weights_dict)
    return(w_gaussian)



def gaussian_moran(scores_df, epsilon = False, threshold = False, visualize = False):
    """
    Calculate Moran's i using Gaussian weight matrix for SPA scores to determine pathway score spatial clustering
    
    Args:
        scores_df (pd.DataFrame): dataframe where columns are pathway names and indices are pixel location
        epsilon (int): epsilon value for gaussian weight matrix. if not given the standard deviation of the pairwise distance matrix is used
        threshold (int): weight threshold below which weight is negligible. (increases speed of computation)
        visualize (bool): If true shows plot of gaussian weight matrix
    Returns:
        Morans_scores (dict): dict with pathway names (key) and moran's i score (value)
    """
    # splitting index
    scores_df.index =  scores_df.index.str.split("_", expand=True)
    # extract x y structure
    sel = pd.DataFrame(scores_df.iloc[:,0])
    coords = list(sel.index)
    #process coords to only include values
    for i, pix in enumerate(coords):
        proc = (int(pix[0][1:]), int(pix[1][1:]))
        coords[i] = proc
    
    coords_array = np.array(coords)

    # Compute full pairwise distance matrix
    distances = cdist(coords_array, coords_array, metric='euclidean')

    # Set a bandwidth parameter for the Gaussian kernel
    # You can choose this manually or based on data spread
    if epsilon == False:
        epsilon = np.std(distances) 

    # Compute Gaussian kernel weights
    w_matrix_full = np.exp(-(distances ** 2) / (2 * epsilon ** 2))

    # Optional: Set diagonal to 0 (no self-weight if you prefer)
    np.fill_diagonal(w_matrix_full, 0)
    # Optional: Row-standardize if needed
    w_matrix_full = w_matrix_full / w_matrix_full.sum(axis=1, keepdims=True)

    # Now w_matrix_full[i, j] gives you the weight between point i and j
    # join 2 levels of index

    scores_df.index = scores_df.index.map(lambda x: f"{x[0]}_{x[1]}")
    if visualize:
        vis = pd.DataFrame(w_matrix_full[200], index=scores_df.index)
        vis = img_conversion(vis)
        plt.imshow(vis)
        plt.title("Gaussian Weight matrix")
        plt.colorbar()

    # convert np array weights to pysal object
    w_gaussian = pysal_w_conv(w_matrix_full, threshold)
    moran_scores = {}
    for pathway in scores_df.columns:
        selected = pd.DataFrame(scores_df[pathway])
        moran = Moran(selected, w_gaussian)
        moran_scores[pathway] = moran.I  # Moran’s I value
    return moran_scores
