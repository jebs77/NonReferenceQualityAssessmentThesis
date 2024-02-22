import csv

# Path to the input CSV file
input_csv_path = 'output_counts.csv'

# Path to the output CSV file
output_csv_path = 'output_sorted.csv'

# Read data from the input CSV file
with open(input_csv_path, 'r', newline='') as infile:
    csvreader = csv.reader(infile)
    data = list(csvreader)

# Sort the data based on the values in the first column
sorted_data = sorted(data, key=lambda x: x[0])

# Write the sorted data to the output CSV file
with open(output_csv_path, 'w', newline='') as outfile:
    csvwriter = csv.writer(outfile)
    csvwriter.writerows(sorted_data)

print(f"CSV file has been sorted alphabetically and written to {output_csv_path}")