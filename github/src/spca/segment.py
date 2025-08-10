import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import cluster


# Visualize PCA with RGB and individually
def pca_rgb_visualize(PCA_df, pc_plot = True):
    """
    Visualizes 3 PC PCA results as an RGB image and optionally each PC individually. Results must come as a pd df with PCs as columns and pixel positions as index
    Args:
        PCA_df (pd.DataFrame): of PCA with PCs as cols
        pc_plot (bool): if True displays each PC individually as well
    Returns:
        pca_rgb:    np array containing max min normalized PCA results in xy format. shape = (y, x, pcs)
    """
    pca_rgb = []
    # convert 3 pcs into a np array of shape (height, width) and append to each other
    for pc in PCA_df.columns:
        pc_series = PCA_df[pc]
        # make vector to 2d matrix by separating x and y tuples
        pc_series.index = pc_series.index.str.split("_", expand = True)   # slick one liner to  turn vector to 2d matrix
        selected_pc = pc_series.reset_index()
        # change coords to int
        selected_pc['level_0'] = selected_pc['level_0'].str[1:].astype(int)
        selected_pc['level_1'] = selected_pc['level_1'].str[1:].astype(int)
        matrix_pc = selected_pc.pivot(index='level_0', columns='level_1', values=pc)
        # max min normalize (0-1)
        norm_matrix_pc = (matrix_pc - matrix_pc.min()) / (matrix_pc.max() - matrix_pc.min())
        # visualize individual PCs.
        if pc_plot == True:
            plt.figure()
            plt.title(str(pc))
            plt.imshow(norm_matrix_pc)
            plt.show()
            plt.close
        pca_rgb.append(np.array(norm_matrix_pc))
    pca_rgb = np.array(pca_rgb)
    pca_rgb = pca_rgb.transpose((1,2,0))
    # visualize as rgb image
    plt.title("3 Principal components visualized as RGB image")
    plt.imshow(pca_rgb)
    return(pca_rgb)


def DBScan_array(array, epsilon, minimun_samples, visualize = True, normalize_coords=True):
    """
    clusters an array of values using DBScan to return a mask with NaN values
    
    Args:
        array:              numpy Array of either PC values, TIC values, or ratio of PC1 loading intensities. Shape (x,y,1 or 3 for PCs)
        epsilon:            The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function. \n
        minimum_samples:    The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself. If min_samples is set to a higher value, DBSCAN will find denser clusters, whereas if it is set to a lower value, the found clusters will be more sparse.\n
        visualize:          Bool which if True visualize the mask
        normalize_coords:   Bool which if True normalizes the input and coordinates between 0 and 1
    Output:
        Mask: np.array same shape as flattened matrix with np.nan values for non-sample tissue
    """

    # Assume selected_matrix is your pivoted DataFrame (level_0, level_1 → TIC values)
    coords = np.array([(x, y) for x in range(array.shape[0]) for y in range(array.shape[1])])  # Pixel coordinates


    # check if 3 PCs, log ratio of loading intensities, or TIC (shape)
    if len(array.shape) > 2:
        ar_shape = array.shape[:-1]
        # Step 1: Reshape the array to combine height and width into a single dimension
        reshaped_array = array.reshape(-1, 3)
        ar_shape = reshaped_array.shape
        # Step 2: Convert the reshaped array to a list of tuples
        values = [tuple(vals) for vals in reshaped_array]

    # if TIC, or log ratio of loading intensities just flatten
    else:
        ar_shape = array.shape
        values = array.flatten().reshape(-1, 1)  # Flatten to 1D
    
    # Stack coordinates & intensities for clustering
    data = np.hstack((coords, values))
    # Normalize intensity values for better clustering
    if normalize_coords:
        data = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0))

    dbscan = cluster.DBSCAN(eps=epsilon, min_samples=minimun_samples)  # Adjust eps & min_samples for best results
    labels = dbscan.fit_predict(data)
    # Convert non_sample points to NaN
    lab_reshape = labels.reshape(ar_shape)
    # # need to get to massage to be image shape again
    # Convert noise points to NaN
    segmented_matrix = np.where(labels.reshape(lab_reshape.shape) < 0, array, np.nan)
    #return segmented_matrix
    if visualize:
        plt.figure()
        plt.title("DBScan segmentation")
        plt.imshow(segmented_matrix)
        plt.close
    return segmented_matrix
