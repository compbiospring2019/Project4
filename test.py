# functions used to create the model
import utils
import os
from random import sample
import sys
from math import exp

empty_row = [-1] * 20
acids_list = ['A', 'C', 'E', 'D', 'G', 'I', 'H', 'K', 'F', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)

# Set the directory for testing output
rr_output_directory = this_script + "/testing-rr-output"

def test(pssm_files, pssm_dir, rr_dir):
    """
    Test the model. Generate contact probabilities and save them in .rr format sorted in descending order.
    """
    model = utils.load_model()
    rr_rows = []
    # create directory for .rr files created by testing
    try:
        os.mkdir(rr_output_directory)
    except FileExistsError:
        print("Removing old .rr files...")
        old_rr_files = utils.read_directory_contents(rr_output_directory, '.rr')
        for old_rr_file in old_rr_files:
            os.remove(os.path.join(rr_output_directory, old_rr_file))
    for pssm_file in pssm_files:
        # get the feature values of every pair
        test_matrix = build_test_matrix(pssm_file, pssm_dir)
        # get the sequence associated with this pssm
        rr = utils.read_rr(pssm_file.replace('.pssm', '.rr'), rr_dir)
        for pair in range(len(test_matrix)):
            contact_probability = calculate_contact_probability(model, test_matrix[pair])
            # create row
            rr_row = [pair[0], pair[1], 0, 8, contact_probability]
            # rr_rows.append(" ".join([str(x) for x in rr_row]))
            rr_rows.append(rr_row)
        # sort rows
        rows.sort(key = lambda row: row[4], reverse = True)
        # write to file
        file_name = os.path.join(rr_output_directory, pssm_file.replace('.pssm', '.rr'))
        with open(file_name, 'w') as file:
            # stringify all rows and put in file
            for row in rr_rows:
                file.write(" ".join([str(x) for x in row]))

def calculate_contact_probability(model, pair):
    """
    :return: the probability that the given pair is in contact
    """
    n = exp(model[0] + sum([model[i + 1] * pair[i] for i in range(200)]))
    return n / (1 + n)

def build_test_matrix(pssm_file, pssm_dir):
    """
    :return: a list of dictionaries with key = (i, j) and value = a list of feature values
    """
    pssm = utils.read_pssm(pssm_file, pssm_dir)
    test_matrix = []
    for i in range(len(pssm)):
        for j in range(i + 5, len(pssm)):
            feature = {}
            row = get_five(pssm, i)
            row.extend(get_five(pssm, j))
            feature[(i, j)] = row
            test_matrix.append(feature)
    return test_matrix

def get_five(pssm, i):
    """
    :return: a list of the 100 feature values in the sliding window around i
    """
    values = []
    for row_offset in range(-2, 3):
        if i + row_offset < 0:
            values.extend(empty_row)
        elif i + row_offset >= len(pssm):
            values.extend(empty_row)
        else:
            values.extend([pssm[i + row_offset][k] for k in acids_list])
    return values

