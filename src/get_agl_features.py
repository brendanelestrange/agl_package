#!/usr/bin/env python

"""
Introduction:
    MWCGS: Multiscale Weighted Colored Subgraphs (C++ Accelerated)

Author:
    Masud Rana (masud.rana@uky.edu)
    (Updated for C++ Integration)
"""

import sys
import os
import pandas as pd
import ntpath
import argparse
import time

import agl_cpp

PROTEIN_TYPES = [
    "C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3", 
    "CG", "CG1", "CG2", "CH2", "CZ", "CZ2", "CZ3", 
    "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", 
    "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OXT", 
    "SD", "SG"
]

LIGAND_ELEMS = [
    "As", "B", "Be", "Br", "C", "Cl", "Co", "Cu", "F", "Fe", 
    "H", "Hg", "I", "Ir", "Mg", "N", "O", "Os", "P", "Pt", 
    "Re", "Rh", "Ru", "S", "Sb", "Se", "Si", "Te", "V", "Zn"
]

# The 10 stats calculated in agl_binding.cpp
STATS = [
    "counts", "sum", "min", "max", "mean", 
    "median", "std_dev", "variance", "size", "sq_norm"
]

class AlgebraicGraphLearningFeatures:

    # Attempt to load kernels, or handle gracefully if missing
    try:
        df_kernels = pd.read_csv('../../utils/kernels.csv')
    except FileNotFoundError:
        df_kernels = pd.DataFrame()

    def __init__(self, args):
        self.kernel_index = args.kernel_index
        self.cutoff = args.cutoff
        self.path_to_csv = args.path_to_csv
        self.matrix_type = args.matrix_type
        self.data_folder = args.data_folder
        self.feature_folder = args.feature_folder
        
        # Generate the column headers once
        self.column_headers = self._generate_headers()

    def _generate_headers(self):
        """
        Replicates the nested loops in C++ to create matching CSV headers.
        Outer Loop: Protein Types
        Inner Loop: Ligand Elements
        Innermost: 10 Stats
        """
        headers = []
        for p_type in PROTEIN_TYPES:
            for l_elem in LIGAND_ELEMS:
                pair_name = f"{p_type}_{l_elem}"
                for stat in STATS:
                    headers.append(f"{pair_name}_{stat}")
        return headers

    def get_agl_features(self, parameters):

        df_pdbids = pd.read_csv(self.path_to_csv)
        pdbids = df_pdbids['PDBID'].tolist()
        pks = df_pdbids['pK'].tolist()

        df_features = pd.DataFrame(columns=self.column_headers)

        print(f"Processing {len(pdbids)} complexes using agl_cpp...")

        for index, _pdbid in enumerate(pdbids):
            lig_file = f'{self.data_folder}/{_pdbid}/{_pdbid}_ligand.mol2'
            pro_file = f'{self.data_folder}/{_pdbid}/{_pdbid}_protein.pdb'

            try:
                feature_values = agl_cpp.get_agl_scores(
                    pro_file,
                    lig_file,
                    float(parameters['cutoff']),
                    float(parameters['power']),
                    float(parameters['tau']),
                    str(parameters['type'])
                )
                
                # Assign to dataframe
                df_features.loc[index] = feature_values

            except Exception as e:
                print(f"Error processing {_pdbid}: {e}")
                # Fill with zeros if C++ fails (e.g., file not found)
                df_features.loc[index] = [0.0] * len(self.column_headers)

        df_features.insert(0, 'PDBID', pdbids)
        df_features.insert(1, 'pK', pks)

        return df_features

    def main(self):
        if not self.df_kernels.empty and self.kernel_index is not None:
            k_row = self.df_kernels.loc[self.kernel_index]
            parameters = {
                'type': k_row['type'],
                'power': k_row['power'],
                'tau': k_row['tau'],
                'cutoff': self.cutoff
            }
        else:
            parameters = {
                'type': 'exponential',
                'power': 6.0,
                'tau': 4.0,
                'cutoff': self.cutoff
            }

        df_features = self.get_agl_features(parameters)

        csv_file_name_only = ntpath.basename(self.path_to_csv).split('.')[0]
        
        k_idx_str = self.kernel_index if self.kernel_index is not None else "custom"
        output_file_name = f'{csv_file_name_only}_agl_{self.matrix_type}_matrix_ker{k_idx_str}_cutoff{self.cutoff}.csv'

        if not os.path.exists(self.feature_folder):
            os.makedirs(self.feature_folder)

        output_path = f'{self.feature_folder}/{output_file_name}'
        df_features.to_csv(output_path, index=False, float_format='%.5f')
        print(f"Saved features to: {output_path}")


def get_args(args):
    parser = argparse.ArgumentParser(description="Get AGL EAT Features")

    parser.add_argument('-k', '--kernel-index', help='Kernel Index (see kernels/kernels.csv)',
                        type=int)
    parser.add_argument('-c', '--cutoff', help='distance cutoff to define binding site',
                        type=float, default=12.0)
    parser.add_argument('-f', '--path_to_csv',
                        help='path to CSV file containing PDBIDs and pK values')
    parser.add_argument('-m', '--matrix_type', type=str,
                        help="type of graph matrix: either 'Laplacian' or 'Adjacency'",
                        default='Laplacian')
    parser.add_argument('-dd', '--data_folder', type=str,
                        help='path to data folder directory')
    parser.add_argument('-fd', '--feature_folder', type=str,
                        help='path to the directory where features will be saved')

    args = parser.parse_args()

    return args


def cli_main():
    args = get_args(sys.argv[1:])

    print(args)
    AGL_Features = AlgebraicGraphLearningFeatures(args)

    AGL_Features.main()


if __name__ == "__main__":

    t0 = time.time()

    cli_main()

    print('Done!')
    print('Elapsed time: ', time.time()-t0)