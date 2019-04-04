# functions used to train and create the model
import utils
import os
from random import sample
import sys
from math import exp

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}
acids_list = ['A', 'C', 'E', 'D', 'G', 'I', 'H', 'K', 'F', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)

# Global variables
SAMPLE_SIZE = 5
STEP_SIZE = 0.001


def main():
    train()


def build_feature_matrix(pssm_files, pssm_dir, rr_dir):
    """
    Builds a feature matrix based on PSSM and RR files
    :return: List of dictionaries. Keys are feature numbers(0 - 199) and 'class'
    """
    print('Building feature matrix...')
    feature_matrix = []

    for pssm_file in pssm_files:
        # For each training file, read in the PSSM matrix and the RR file
        pssm = utils.read_pssm(pssm_file, pssm_dir)
        rr = utils.read_rr(pssm_file.replace('.pssm', '.rr'), rr_dir)
        intermediate_matrix = build_small_matrix(pssm)

        # Get a random sample of (i, j) pairs to balance data
        pairs = [(i, j) for i in range(len(intermediate_matrix)) for j in range(i+5, len(intermediate_matrix)) if (i, j) not in rr]
        pairs = sample(pairs, len(intermediate_matrix))

        def build_row(class_label):
            # Build row for feature matrix:
            feature_row = {'class': class_label}
            feature_row.update(intermediate_matrix[i])
            for feat_num in range(len(intermediate_matrix[j])):
                feature_row[100 + feat_num] = intermediate_matrix[j][feat_num]
            # feature_matrix.append(feature_row)
            return feature_row

        # For each i, j pair, make a row for the feature matrix
        for i, j in pairs:
            feature_matrix.append(build_row(0))
        for i, j in rr.keys():
            feature_matrix.append(build_row(1))
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


def train():
    """
    Train the model using gradient ascent. Save the model.
    """
    # Build feature matrix
    pssm_list, rr_list, pssm_dir, rr_dir = parse_args()
    feature_matrix = build_feature_matrix(pssm_list, pssm_dir, rr_dir)

    w_vector = new_w_vector()
    gradient_vector = None

    print('Training the model...')
    while not reached_top(w_vector, gradient_vector):
        gradient_vector = calc_gradient(w_vector, feature_matrix)
        w_vector = update_w(w_vector, gradient_vector)

    # Save the model to the file
    utils.write_model(w_vector)


def calc_gradient(w_vector, matrix):
    """
    Calculates the gradient based on (SAMPLE_SIZE) training examples
    :return: gradient vector
    """
    gradient_vector = [0] * len(w_vector)
    training_data = sample(matrix, SAMPLE_SIZE)

    for training_example in training_data:
        # Calculate P(Y=1|X,w)
        sum_of_w = w_vector[0]
        for i in range(1, len(w_vector)):
            sum_of_w += w_vector[i] * training_example[i-1]
        p_hat = 1.0 / (1 + exp(sum_of_w))

        # Deal with w0
        gradient_vector[0] += training_example['class'] - p_hat

        # For each feature, calculate the gradient
        for i in range(1, len(w_vector)):
            gradient_vector[i] += training_example[i-1] * (training_example['class'] - p_hat)

    return gradient_vector


def update_w(w_vector, gradient_vector):
    """
    Updates each w value in the w_vector
    :return: w_vector
    """
    for index in range(len(w_vector)):
        w_vector[index] += STEP_SIZE * gradient_vector[index]
    return w_vector


def reached_top(w_vector, gradient_vector):
    """
    Check if we've reached the top of the mountain
    :return: boolean
    """
    if not gradient_vector:
        return False
    return True


def new_w_vector():
    """
    Make a new w vector.
    TODO: Random values for 'base camp'?
    :return: list of length 201
    """
    return [1] * 201


err_msg = '''
Please enter two directory names (absolute paths)
containing sequences for linear regression training data
(with double quotes around them if they have spaces).
The directory with PSSM files should come first, 
followed by the path to the .rr files.'''


def parse_args():
    if len(sys.argv) < 3:
        print(err_msg)
        sys.exit()

    try:
        # Get the lists of pssm and rr file names
        pssm = utils.read_directory_contents(sys.argv[1], '.pssm')
        rr = utils.read_directory_contents(sys.argv[2], '.rr')
    except:
        # Given paths are not valid directories
        print(err_msg)
        sys.exit()

    # Return list of pssm & rr files, and their parent directories
    return pssm, rr, sys.argv[1], sys.argv[2]


if __name__ == '__main__':
    main()
