import pandas as pd
import matplotlib.pyplot as plt
import os



col_names = ["percent_reads", "reads_rooted", "reads_taxon", "rank", "tax_id", "tax_name"]

df = pd.read_csv('all_samples_genus.tsv', sep='\t', names=col_names)

bacterial_genera = pd.read_csv('bacterial_genera.tsv', sep='\t', header=0)['genus'].tolist()

df = df.groupby(['tax_id','tax_name']).agg({'reads_rooted':'sum'}).reset_index()


df = df[df['tax_name'].isin(bacterial_genera)].sort_values('reads_rooted', ascending=False)

df.to_csv('genera_ranking.tsv', sep='\t', index=False, header=False)


