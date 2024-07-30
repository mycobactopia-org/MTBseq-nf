import sys
import csv

if len(sys.argv) < 2:
    print("Uso: python make_symmetric_matrix.py <nome_do_arquivo>")
    sys.exit(1)


input_file = sys.argv[1]

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
        if i==j:
            matrix[i + 1][j + 1] = 0
        if value: 
            matrix[i + 1][j + 1] = value
            matrix[j + 1][i + 1] = value 
        

with open('symmetric_matrix.txt', 'w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(matrix)