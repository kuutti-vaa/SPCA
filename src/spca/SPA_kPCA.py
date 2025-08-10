import pandas as pd
import numpy as np
from sklearn.decomposition import KernelPCA
import sspa.utils as utils
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.estimator_checks import check_estimator
from sklearn.base import BaseEstimator

class spa_KPCA(BaseEstimator):
    """
    Kernel PCA method for spatial pathway analysis (SPA)

    Args:
        pathway_df: (pd.DataFrame): pandas DataFrame of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        formula_chebis (dict): dictionary of dataset formulas (keys) and a list of corresponding possible ChEBI ids (values)
        min_entity (int): minimum number of metabolites mapping to pathways for SPA to be performed
    """
    def __init__(self, pathway_df, formula_chebis, min_entity=2, random_state = 0):
        self.pathway_df = pathway_df 
        self.pathway_dict = utils.pathwaydf_to_dict(pathway_df)
        self.formula_chebis = formula_chebis
        self.min_entity = min_entity
        self.pathways_filt = {}
        self.pathway_matrices = []
        self.fitted_models = []
        self.pathway_ids = []
        self.pathway_names = []
        self.pathway_formulas = []
        self.random_state = random_state

    def fit(self, X, y = None):
        """
        Fit the model with X
        Args:
            X (pd.DataFrame): pandas DataFrame containing MSI data containing pixels (rows) and entities (columns)
            Do not include metadata columns
            Returns: 
            self : object
        """
        self.X_ = X
        self.y_ = y
        for pathway_id, compounds in self.pathway_dict.items():
            try:
                pathway_metabs = set(self.pathway_dict[pathway_id])

                # # Use list comprehension with 'any' to check for intersections without converting each v to a set
                pathway_metabs = [k for k, v in self.formula_chebis.items() if any(chebi in pathway_metabs for chebi in v)]

            except KeyError:
                continue
            if not pathway_metabs or len(pathway_metabs) < self.min_entity:
                continue
            else:
                self.pathway_ids.append(pathway_id)
                self.pathway_names.append(self.pathway_df.loc[pathway_id, 'Pathway_name'])
                self.pathway_formulas.append(pathway_metabs)
                pathway_mat = X.loc[:, pathway_metabs]
                self.pathway_matrices.append(pathway_mat)
                kpca = KernelPCA(n_components=2, kernel="rbf", random_state=self.random_state)
                self.fitted_models.append(kpca.fit(pathway_mat.to_numpy()))

                print(f"pathway {self.pathway_df.loc[pathway_id, 'Pathway_name']} SPA complete with {len(pathway_metabs)} metabs")
        self.pathways_filt = {k: v for k, v in self.pathway_df.items() if k in self.pathway_ids}
        self.is_fitted_ = True
        return self
    
    def transform(self, X, y = None):
        """
                Transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame containing MSI data containing pixels (rows) and entities (columns)
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the kPCA method. Columns represent pathways and rows represent samples.
        """
        check_is_fitted(self, "is_fitted_")
        scores = []
        #  here massage data into clean pd.df
        for i,model in enumerate(self.fitted_models):
            mat = self.pathway_matrices[i]
            data = model.transform(mat.to_numpy())
            scores.append(data[:,0])
        scores_df = pd.DataFrame(np.array(scores).T, columns=self.pathway_names, index = X.index)
        return scores_df
    

    def fit_transform(self, X, y=None):
        """
        Fit the model with X and transform X.

        Args:
            X (pd.DataFrame): pandas DataFrame containing MSI data containing pixels (rows) and entities (columns)
            Do not include metadata columns
            Returns: 
            pandas DataFrame of pathway scores derived using the kPCA method. Columns represent pathways and rows represent samples.
        """
        self.X_ = X
        self.y_ = y

        self.fit(X)
        return self.transform(X)