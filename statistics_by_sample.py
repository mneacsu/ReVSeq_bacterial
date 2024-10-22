import pandas as pd
import os

dir_path = 'kraken2_results'

top_20_genera = pd.read_csv('genera_ranking.tsv', sep='\t', header=None).iloc[:20,1].tolist()
pca_df = pd.DataFrame(columns=top_20_genera)

col_names = ["percent_reads", "reads_rooted", "reads_taxon", "rank", "tax_id", "tax_name"]

bacterial_genera = pd.read_csv('bacterial_genera.tsv', sep='\t', header=0)['genus'].tolist()

pca_df = pd.DataFrame(columns=top_20_genera)

for plate in os.listdir(dir_path):
	for sample in os.listdir(os.path.join(dir_path, plate)):

		df = pd.read_csv(os.path.join(dir_path, plate, sample, sample + '.tsv'), sep='\t', names = col_names)

		bacterial_reads = df[df['tax_name']=='Bacteria']['reads_rooted'].tolist()[0] if not df[df['tax_name']=='Bacteria'].empty else 0
		
		df['percent_bacteria'] = df['reads_rooted']/bacterial_reads if bacterial_reads else 0

		df = df[df['tax_name'].isin(bacterial_genera)]

		df= df.set_index('tax_name')

		perc_per_genus = [df.loc[genus,'percent_bacteria'] if (genus in df.index) else 0 for genus in top_20_genera]
        
		pca_df.loc[sample,:] = perc_per_genus


pca_df.to_csv('genera_by_sample.csv')