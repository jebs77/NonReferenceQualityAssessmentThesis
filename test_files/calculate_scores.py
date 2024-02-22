import csv
import os

directory = 'VSENSE\scores'

words_to_filter = ['longdress', 'loot', 'soldier','redandblack', 'orig']

column_counts = {}

# Iterate through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        filepath = os.path.join(directory, filename)
        
        # Open and read the CSV file
        with open(filepath, 'r', newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            
            for row in csvreader:
                if row:  
                    first_column_value = row[0]
                    if first_column_value not in words_to_filter:  
                        if first_column_value in column_counts:
                            column_counts[first_column_value] += 1
                        else:
                            column_counts[first_column_value] = 1

# Write the counts to a CSV file
output_file_path = 'output_counts_filtered.csv'
with open(output_file_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Value', 'Count'])
    for value, count in column_counts.items():
        csvwriter.writerow([value, count])

print(f"Filtered counts have been written to {output_file_path}")