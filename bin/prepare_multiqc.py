#!/usr/bin/env python3

import sys
import csv
import argparse

def process_matrix_file(input_file, output_file):
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        data = list(reader)

    sample_names = [row[0] for row in data]

    n = len(data)
    matrix = [['' for _ in range(n + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        matrix[i][0] = sample_names[i - 1]
        matrix[0][i] = sample_names[i - 1]

    for i, row in enumerate(data):
        for j, value in enumerate(row[1:]):
            if i == j:
                matrix[i + 1][j + 1] = 0
            if value:
                matrix[i + 1][j + 1] = value
                matrix[j + 1][i + 1] = value

    with open(output_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(matrix)

def process_tab_file(input_file, output_file):
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        data = [row[1:] for row in reader]

    data = [[value.replace("'", "") for value in row] for row in data]

    with open(output_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(data)

def process_groups_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    start_index = 0
    for i, line in enumerate(lines):
        if "### Output as lists:" in line:
            start_index = i + 1
            break

    lines = lines[start_index:]
    lines.insert(0, "Sample\tClustering\tGroup\n")

    with open(output_file, 'w') as file:
        file.writelines(lines)

def main():
    parser = argparse.ArgumentParser(description='Prepares tab files for multiqc importation with tsv.')
    parser.add_argument('--mapping', required=True, help='Mapping_and_Variant_Statistics.tab')
    parser.add_argument('--strain', required=True, help='Strain_Classification.tab')
    parser.add_argument('--groups', required=True, help='Groups file')
    parser.add_argument('--matrix', required=True, help='.matrix')

    args = parser.parse_args()

    process_tab_file(args.mapping, 'Mapping_and_Variant_Statistics_processed.tsv')
    process_tab_file(args.strain, 'Strain_Classification_processed.tsv')
    process_groups_file(args.groups, 'Groups_processed.tsv')
    process_matrix_file(args.matrix, 'symmetric_matrix.txt')

if __name__ == "__main__":
    main()
