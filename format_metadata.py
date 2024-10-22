import pandas as pd
import numpy as np

df = pd.read_csv('mara_all_data.csv', header=0) #sep=';', encoding='latin-1'
df['strains_PCR'] = df['strains_PCR'].astype(str)

df_fmt = pd.DataFrame(columns = ['sample', 'substrains', 'strains_seq', 'strains_PCR', 'multiple_seq', 'multiple_PCR', 'positive_seq', 'positive_PCR', 'perfect_match', 'partial_match', 'plate', 'aligned_reads', 'prescriber', 'ent_date', 'panels_PCR', 'coverage', 'DP5', 'DP10', 'DP20'])


for sample in df['pseudonymized_id'].unique():

	df_sample = df[df['pseudonymized_id']==sample].sort_values('DP10', ascending=False)
	substrains = df_sample[df_sample['DP10']>.2]['substrain_name'].tolist()
	strains_seq = df_sample[df_sample['DP10']>.2]['strain_name'].unique().tolist()
	strains_PCR = df_sample['strains_PCR'].iloc[0].split(',') if df_sample['strains_PCR'].iloc[0] != 'nan' else []
	positive_seq = (len(strains_seq) > 0)
	positive_PCR = (len(strains_PCR) > 0)
	multiple_seq = (len(strains_seq) > 1)
	multiple_PCR = (len(strains_PCR) > 1)
	perfect_match = (set(strains_seq) == set(strains_PCR))
	partial_match = (len(set(strains_seq).intersection(strains_PCR)) > 0)
	plate = df_sample['barcode'].iloc[0]
	aligned_reads = max(df_sample['aligned'])
	prescriber = df_sample['prescriber'].iloc[0]
	ent_date = df_sample['ent_date'].iloc[0]
	panels_PCR = df_sample['panels_PCR'].iloc[0]
	coverage = max(df_sample['coverage'])
	DP5 = max(df_sample['DP5'])
	DP10 = max(df_sample['DP10'])
	DP20 = max(df_sample['DP20'])

	df_fmt = df_fmt._append({'sample': sample, 'substrains': ','.join(substrains), 'strains_seq': ','.join(strains_seq), 'strains_PCR': ','.join(strains_PCR), 'multiple_seq': multiple_seq, 'multiple_PCR': multiple_PCR, 'positive_seq': positive_seq, \
		'positive_PCR': positive_PCR, 'perfect_match': perfect_match, 'partial_match': partial_match, 'plate': plate, 'aligned_reads': aligned_reads, 'prescriber': prescriber, 'ent_date': ent_date, 'panels_PCR': panels_PCR, 'coverage': coverage, 'DP5': DP5, 'DP10': DP10, 'DP20': DP20}, ignore_index = True)

df_fmt.to_csv('all_plates_fmt.csv', index=False)
