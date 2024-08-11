#!/usr/bin/env python3

import sys
import os
import csv
import argparse


def find_file_in_subfolders(filename, root_folder='.'):
    for root, dirs, files in os.walk(root_folder):
        if filename in files:
            return os.path.join(root, filename)
    return None


def process_matrix_file(input_file, output_file):

    # Find the input file in all subfolders
    input_file_path = find_file_in_subfolders(input_file)
    if not input_file_path:
        print(f"File {input_file} not found in any subfolder.")
        sys.exit(1)


    # Read the TSV file into a list of lists
    with open(input_file_path, 'r') as file:
        data = file.readlines()

    # Extract row labels and matrix values
    row_labels = []
    matrix_values = []

    for row in data:
        values = row.strip().split('\t')
        row_labels.append(values[0])
        matrix_values.append(values[1:])

    # Initialize the matrix
    n = len(matrix_values)  # Number of rows/columns
    matrix = [['' for _ in range(n)] for _ in range(n)]

    # Fill the matrix
    for i, row in enumerate(matrix_values):
        for j, value in enumerate(row):
            if value:  # Skip empty values
                matrix[i][j] = value
                matrix[j][i] = value  # Ensure symmetry

    # Write the symmetric matrix to a new TSV file
    with open(output_file, 'w') as file:
        file.write('\t' + '\t'.join(row_labels) + '\n')
        for i, row in enumerate(matrix):
            file.write(row_labels[i] + '\t' + '\t'.join(row) + '\n')


def process_tab_file(input_file, output_file):

    # Find the input file in all subfolders
    input_file_path = find_file_in_subfolders(input_file)
    if not input_file_path:
        print(f"File {input_file} not found in any subfolder.")
        sys.exit(1)


    # Read the TSV file into a list of lists
    with open(input_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        data = [row[1:] for row in reader]

    data = [[value.replace("'", "") for value in row] for row in data]

    with open(output_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(data)

def process_groups_file(input_file, output_file):


    # Find the input file in all subfolders
    input_file_path = find_file_in_subfolders(input_file)
    if not input_file_path:
        print(f"File {input_file} not found in any subfolder.")
        sys.exit(1)


    # Read the TSV file into a list of lists
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    start_index = 0
    for i, line in enumerate(lines):
        if "### Output as lists:" in line:
            start_index = i + 1
            break

    lines = lines[start_index:]
    lines.insert(0, "Sample\tGroup\n")

    with open(output_file, 'w') as file:
        file.writelines(lines)

def main():
    parser = argparse.ArgumentParser(description='Prepares tab files for MultiQC.')
    parser.add_argument('--mapping', required=True, help='Mapping_and_Variant_Statistics file')
    parser.add_argument('--strain', required=True, help='Strain_Classification file')
    parser.add_argument('--groups', required=True, help='Groups file')
    parser.add_argument('--matrix', required=True, help='Groups matrix file')

    args = parser.parse_args()

    process_tab_file(args.mapping, 'Mapping_and_Variant_Statistics_processed.tsv')
    process_tab_file(args.strain, 'Strain_Classification_processed.tsv')
    process_groups_file(args.groups, 'Groups_processed.tsv')
    process_matrix_file(args.matrix, 'Matrix_processed.tsv')

if __name__ == "__main__":
    main()
