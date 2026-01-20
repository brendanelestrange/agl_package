import csv

input_file = "PDBbind_core_set_v2007.2.lst"
output_file = "CASF_2007_CoreSet.csv"

data = []

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue  # skip comments and blank lines
        parts = line.split()
        if len(parts) < 4:
            continue  # skip malformed lines
        pdbid = parts[0]
        pk = parts[3]
        data.append([pdbid, float(pk)])

# Write to CSV
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["index", "PDBID", "pK"])
    for i, row in enumerate(data):
        writer.writerow([i] + row)

print(f"Converted {len(data)} entries to {output_file}")
