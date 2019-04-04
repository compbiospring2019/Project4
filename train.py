# functions used to train and create the model
import utils
import os

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}
acids_list = ['A', 'C', 'E', 'D', 'G', 'I', 'H', 'K', 'F', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


def build_feature_matrix(pssm_files, pssm_dir, rr_dir):
    """
    Builds a feature matrix based on PSSM and RR files
    :return: List of dictionaries. Keys are feature numbers(0 - 199) and 'class'
    """
    feature_matrix = []

    for pssm_file in pssm_files:
        # For each training file, read in the PSSM matrix and the RR file
        pssm = utils.read_pssm(pssm_file, pssm_dir)
        rr = utils.read_rr(pssm_file.replace('.pssm', '.rr'), rr_dir)
        intermediate_matrix = build_small_matrix(pssm)

        # For each i, j pair, make a row for the feature matrix
        for i in range(len(intermediate_matrix)):
            for j in range(i+5, len(intermediate_matrix)):
                # Build row for feature matrix:
                feature_row = {'class': 0}
                feature_row.update(intermediate_matrix[i])
                for feat_num in range(len(intermediate_matrix[j])):
                    feature_row[100 + feat_num] = intermediate_matrix[j][feat_num]

                # If (i, j) is in the .rr file, they're in contact
                if (i, j) in rr:
                    feature_row['class'] = 1
                feature_matrix.append(feature_row)
    return feature_matrix


def build_small_matrix(pssm):
    """
    Builds the 'small' intermediate matrix. Rows are 100 features from sliding window of size 5.
    :return: List of dictionaries. Keys are feature numbers (0 - 99)
    """
    small_matrix = []
    for row_num in range(len(pssm)):
        # For each amino acid in the PSSM, build a line for the feature matrix
        feature = {}
        for row_offset in range(-2, 3):
            if row_num + row_offset < 0:
                # We're at the top of the PSSM
                values = empty_row
            elif row_num + row_offset >= len(pssm):
                # We're at the bottom of the PSSM
                values = empty_row
            else:
                # We're somewhere in the middle
                values = pssm[row_num + row_offset]
            for val_num, acid in enumerate(acids_list):
                feature[((row_offset + 2) * 20) + val_num] = values[acid]
        small_matrix.append(feature)
    return small_matrix
