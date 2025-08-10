import pandas as pd
import numpy as np
from sklearn.decomposition import KernelPCA
from sklearn.utils.validation import check_is_fitted
from sklearn.base import BaseEstimator

class sca_KPCA(BaseEstimator):
    """
    Kernel PCA method for spatial class analysis (SCA)

    Args:
        uni_parents: list of unique parent molecules (classes) 
        parent_formula (dict): dict containing parent molecules (classes) and the annotation molecular formulas which map to them
        min_entity (int): minimum number of metabolites mapping to pathways for SCA to be performed
    """
    def __init__(self, uni_parents, parent_formula, min_entity=2, random_state = 0):
        self.uni_parents = uni_parents 
        self.parent_formula = parent_formula
        self.min_entity = min_entity
        self.class_matrices = []
        self.fitted_models = []
        # self.class_ids = []
        self.class_names = []
        # self.class_formulas = []
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
        for parent_mol in self.uni_parents:
            try:
                parent_metabs = self.parent_formula[parent_mol]

            except KeyError:
                continue
            if not parent_metabs or len(parent_metabs) < self.min_entity:
                continue
            else:
                # self.pathway_ids.append(pathway_id)
                self.class_names.append(parent_mol)
                # self.pathway_formulas.append(pathway_metabs)
                class_mat = X.loc[:, parent_metabs]
                self.class_matrices.append(class_mat)
                kpca = KernelPCA(n_components=2, kernel="rbf", random_state=self.random_state)
                self.fitted_models.append(kpca.fit(class_mat.to_numpy()))

                print(f"parent class {parent_mol} SCA complete with {len(parent_metabs)} metabs")
            
        # self.pathways_filt = {k: v for k, v in self.pathway_df.items() if k in self.pathway_ids}
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
            mat = self.class_matrices[i]
            data = model.transform(mat.to_numpy())
            scores.append(data[:,0])
        scores_df = pd.DataFrame(np.array(scores).T, columns=self.class_names, index = X.index)
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