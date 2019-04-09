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
rr_output_directory = parent_directory + "/testing-rr-output"

def main(pssm_files, pssm_dir, rr_dir):
    test(pssm_files, pssm_dir, rr_dir)
    accuracy(rr_dir)

def test(pssm_files, pssm_dir, rr_dir):
    """
    Test the model. Generate contact probabilities and save them in .rr format sorted in descending order.
    """
    model = utils.read_model()
    rr_rows = []

    # create directory for .rr files created by testing
    try:
        os.mkdir(rr_output_directory)
    except FileExistsError:
        print("Removing old .rr files...")
        old_rr_files = utils.read_directory_contents(rr_output_directory, '.rr')
        for old_rr_file in old_rr_files:
            os.remove(os.path.join(rr_output_directory, old_rr_file))

    # test
    print("Testing...")
    for pssm_file in pssm_files:
        # get the feature values of every pair
        test_matrix = build_test_matrix(pssm_file, pssm_dir)
        # get the sequence associated with this pssm
        sequence = utils.read_rr(pssm_file.replace('.pssm', '.rr'), rr_dir)['sequence']
        for pair in range(len(test_matrix)):
            contact_probability = calculate_contact_probability(model, test_matrix[pair])
            # create row
            i, j = list(test_matrix[pair].keys())[0][0], list(test_matrix[pair].keys())[0][1]
            rr_row = [i, j, 0, 8, contact_probability]
            # rr_rows.append(" ".join([str(x) for x in rr_row]))
            rr_rows.append(rr_row)

        # sort rows
        rr_rows.sort(key = lambda row: row[4], reverse = True)

        # write to file
        file_name = os.path.join(rr_output_directory, pssm_file.replace('.pssm', '.rr'))
        with open(file_name, 'w') as file:
            file.write(sequence + "\n")
            # stringify all rows and put in file
            for row in rr_rows:
                file.write(" ".join([str(x) for x in row]) + "\n")

def accuracy(rr_dir):
    """
    Determine the average accuracy of the L/10, L/5, and L/2 most probable contact predictions.
    """
    # metrics
    l10_acc, l5_acc, l2_acc = 0.0, 0.0, 0.0

    # get the .rr files created in testing (predictions)
    rr_outputs = utils.read_directory_contents(rr_output_directory, '.rr')

    # determine the accuracy of each individual .rr file
    num_files = 0
    for rr_output in rr_outputs:
        num_files += 1
        # create a list of the predictions
        predictions = utils.read_rr(rr_output, rr_output_directory)
        predictions_list = listify_rr(predictions)
        predictions_list.sort(key = lambda row: row[2], reverse = True)
        # create a list of the actual contact pairs
        contact_pairs = utils.read_rr(rr_output, rr_dir)
        contact_pairs_list = listify_rr(contact_pairs)

        # metrics
        l10_denom, l10_num, l5_denom, l5_num, l2_denom, l2_num = 0, 0, 0, 0, 0, 0

        L = len(contact_pairs['sequence'])
        for i in range(int(L / 2)):
            # make prediction
            predicted_contact = predictions_list[i][2] > 0.5
            actual_contact = (predictions_list[i][0], predictions_list[i][1]) in contact_pairs

            # update metrics
            if i <= L / 10 and actual_contact:
                l10_denom += 1
                if predicted_contact:
                    l10_num += 1
            if i <= L / 5 and actual_contact:
                l5_denom += 1
                if predicted_contact:
                    l5_num += 1
            if actual_contact:
                l2_denom += 1
                if predicted_contact:
                    l2_num += 1

        # aggregate
        if l10_denom != 0:
            l10_acc += l10_num / l10_denom
        if l5_denom != 0:
            l5_acc += l5_num / l5_denom
        if l2_denom != 0:
            l2_acc += l2_num / l2_denom

    # average
    l10_acc /= num_files
    l5_acc /= num_files
    l2_acc /= num_files
    print("Accuracy")
    print("--------")
    print("L10: " + str(l10_acc))
    print("L5: " + str(l5_acc))
    print("L2: " + str(l2_acc))

def listify_rr(rr):
    pair_keys = [key for key in rr.keys() if key != 'sequence']
    return [[pair[0], pair[1], rr[pair]] for pair in pair_keys]

def calculate_contact_probability(model, pair):
    """
    :return: the probability that the given pair is in contact
    """
    feature_values = list(pair.values())[0]
    n = exp(model[0] + sum([model[i + 1] * feature_values[i] for i in range(200)]))
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

if __name__ == "__main()__":
    main()
