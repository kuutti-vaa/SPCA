import numpy as np
import pandas as pd


def process_image_set(image_set, results):
    """
    Converts image set list of images into pd.dataframe with index as compound formulas, a column containing chebi ids and flattened pixel intensities.
    
    Args:
    image_sets: list of images (from metaspace isotope_image_sets method)
    results: dataframe containing dataset entity molecular formulas which has been filtered by FDR
    returns:
    pixel_intensities: pd.dataframe with mol formula as index, chebi ids as first columns and the matching pixel intensities as the following columns
    """
    # input should be isotope image sets (principal peak images for each feature with fdr less than 0.1)
    positions = []

    for y,vy in enumerate(image_set[0][0]):        # maybe flip it ;)
        for x, vx in enumerate(vy):
            t_str = "y" + str(y) + "_" + "x" + str(x)
            positions.append(t_str)

    # flatten the principal peak of all features
    intensities = np.array([imgs[0].flatten() for imgs in image_set])
    pos_intensities = pd.DataFrame(intensities, columns=positions)

    # add mol formula
    mol_formula = [str(str(imgs.formula) + imgs.adduct[1:]) for imgs in image_set]
    pos_intensities["mol_formula"] = mol_formula
    pos_intensities = pos_intensities[ ["mol_formula"] + [ col for col in pos_intensities.columns if col != "mol_formula" ] ]

    pos_intensities = pos_intensities[pos_intensities["mol_formula"].isin(results["mol_formula"])]
    # add chebi column to pixel intensities df 
    pixel_intensities = pos_intensities
    pixel_intensities["chebi"] = pixel_intensities.mol_formula.map(dict(zip(results["mol_formula"], results["moleculeIds"].values)))

    # move chebi to first column and set index to molecular formula
    pixel_intensities = pixel_intensities[ ["chebi"] + [ col for col in pixel_intensities.columns if col != "chebi" ] ]
    pixel_intensities = pixel_intensities.set_index('mol_formula')
    return pixel_intensities

def SPCA_processing(pixel_intensities, sample_mask = None):
    """ 
    Massages pixel intensity dataframe and log normalizes the intensities. Massaging includes dropping Chebi id column and transposing.
    
    Args:
    pixel_intensities (pd.DataFrame): Dataframe with molecular formulas as index, ChEBI ids as first row and pixel intensities as subsequent rows
    sample_mask (pd.DataFrame): optional argument, if given will drop non-sample pixels
    """
    pixel_intensities_sspa = pixel_intensities.iloc[:, 1:]
    pixel_intensities_sspa = pixel_intensities_sspa.T         
    pixel_intensities_sspa = pixel_intensities_sspa + 1
    pixel_intensities_sspa = pd.DataFrame(np.log2(pixel_intensities_sspa), 
                                            index=pixel_intensities_sspa.index, 
                                            columns=pixel_intensities_sspa.columns)
    # remove columns (features) with all zeros
    pixel_intensities_sspa = pixel_intensities_sspa.loc[:, (pixel_intensities_sspa != 0).any(axis=0)]
    if sample_mask is not None:
        # drop all "non-tissue" pixels
        pixel_intensities_sspa["mask"] = sample_mask["Mask"]
        pixel_intensities_sspa = pixel_intensities_sspa.dropna()
        pixel_intensities_sspa = pixel_intensities_sspa.drop("mask", axis=1)
    return pixel_intensities_sspa


def img_conversion(selected_pathway):
    """ 
    Convert vector of single pathway pixel pseudo-activities into 2d matrix for visalization
    
    Args:
    selected_pathway:   pandacore series containing x position and y position as a str and pixel intensity
    returns:
    matrix:             2d matrix containing pixel values corresponding to x and y 
    """

    # make vector to 2d matrix by separating x and y tuples
    selected_pathway.index = selected_pathway.index.str.split("_", expand = True)   # slick one liner to  turn vector to 2d matrix
    selected_pathway = selected_pathway.reset_index()

    # change coords to int
    selected_pathway['level_0'] = selected_pathway['level_0'].str[1:].astype(int)
    selected_pathway['level_1'] = selected_pathway['level_1'].str[1:].astype(int)
    selected_matrix = selected_pathway.pivot(index='level_0', columns='level_1')
    return selected_matrix
