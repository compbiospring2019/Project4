# Utils for Project 4
import os
from random import sample
import json

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


# Read a biological sequence or RSA sequence from a file:
def read_sequence(file_path, dir=None):
    if file_path is None:
        return None
    sequence = ''
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        for line in f:
            sequence += line.strip()

    return sequence.upper()


def read_directory_contents(path, file_extension):
    if not os.path.isdir(path):
        # This is not a valid directory!
        raise Exception('Not a valid directory!')

    # Return a list of files with the file extension
    ls_dir = os.listdir(path)
    return [file_name for file_name in ls_dir if file_name.endswith(file_extension)]


def read_pssm(file_path, dir=None):
    if dir:
        file_path = os.path.join(dir, file_path)
    pssm = []
    with open(file_path, 'r') as f:
        # Ignore the title line
        title = f.readline()
        if title in ['', '\n']:
            title = f.readline()

        # Get the list of amino acids on the top
        headers = f.readline().strip().split()[:20]

        # Now, read each line of the matrix into a dictionary
        for line in f:
            if line in ['', '\n']:
                break
            line_list = line.strip().split()[:22]
            row = {'this-acid': line_list[1]}
            for acid_num in range(len(headers)):
                row[headers[acid_num]] = int(line_list[acid_num + 2])
            pssm.append(row)

    # Returns a list of dictionaries, where each dict is a row of the matrix
    return pssm


def read_rr(file_path, dir=None):
    if dir:
        file_path = os.path.join(dir, file_path)
    rr = {}
    with open(file_path, 'r') as f:
        # Read in the sequence from the top of the .rr file
        sequence = ''
        line = f.readline().strip()
        while '' not in line:
            sequence += line
            line = f.readline().strip()

        # rr['sequence'] = sequence

        # Add each line to the dictionary
        for line in f:
            if line in ['', '\n']:
                break
            rr.update(parse_rr_line(line))

    # Returns a dictionary with tuples as keys, and intra residue distance as value
    return rr


def parse_rr_line(rr_line):
    """
    Splits a line of the .rr file into i, j, and distance,
    and returns a dictionary with key = (i, j) and value = distance
    """
    rr_line = rr_line.strip()
    rr_parts = rr_line.split()
    i = int(rr_parts[0])
    j = int(rr_parts[1])
    distance = float(rr_parts[4])
    rr_dict = {(i-1, j-1): distance}
    return rr_dict


# prior - probability of class label
# dists - list of means and standard deviations
def read_dist(file_path, dir=None):
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        prior = float(f.readline())
        dists = [[float(x) for x in line] for line in f]
    return prior, dists


# divide the .pssm files into training and testing sets
def split_files(pssm_list, rr_list):
    test_correct_pssm_files(pssm_list, rr_list)
    pssm_train = sample(pssm_list, int(0.75 * len(pssm_list)))
    pssm_test = [pssm_name for pssm_name in pssm_list if pssm_name not in pssm_train]
    return pssm_train, pssm_test


# make sure each .pssm has a corresponding .rr
def test_correct_pssm_files(pssm_list, rr_list):
    for pssm_name in pssm_list:
        if pssm_name.replace('.pssm', '.rr') not in rr_list:
            raise Exception('PSSM files don\'t match up with .rr files: {}'.format(pssm_name))


def write_model(model, file_name='model.json', dir=parent_directory):
    """
    Write JSON object to a file
    :return: None
    """
    if dir:
        file_name = os.path.join(dir, file_name)
    with open(file_name, 'w') as outfile:
        json.dump(model, outfile)


def read_model(file_name='model.json', dir=parent_directory):
    """
    Reads in the JSON object from model.json
    :return: Loaded contents of model file
    """
    if dir:
        file_name = os.path.join(dir, file_name)
    with open(file_name, 'r') as file:
        model = json.load(file)

    return model
