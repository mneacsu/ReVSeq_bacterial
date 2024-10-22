import pandas as pd
import os 

dir_path = 'kraken2_results'
col_names = ["percent_reads", "reads_rooted", "reads_taxon", "rank", "tax_id", "tax_name"]
genera =[]

for plate in os.listdir(dir_path):
	for sample in os.listdir(os.path.join(dir_path, plate)):

		df = pd.read_csv(os.path.join(dir_path, plate, sample, sample + '.tsv'), sep = '\t', names = col_names)

		df = df[df['rank']=='G']

		genera = genera + df['tax_name'].tolist()

pd.DataFrame({'genus': genera}).to_csv('detected_genera.csv', index=False)