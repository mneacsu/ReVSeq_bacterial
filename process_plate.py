import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("kraken2_dir", help="directory where the kraken2 reports are located")
parser.add_argument("reads_dir", help="directory where the raw reads are located")
parser.add_argument("metadata", help="metadata table")
parser.add_argument("output_dir", help="output directory")
parser.add_argument("--search_VF", help="option to search for virulence factors", action="store_true")
parser.add_argument("-d", "--VF_database", help="path to the virulence factor database", type=str, default="out.vfdb.prepareref")
parser.add_argument("-t", "--threads", help="number of threads", type=int, default=1)


args = parser.parse_args()
    
path_to_kraken2_results = args.kraken2_dir
path_to_raw_reads = args.reads_dir

metadata_file = args.metadata

search_AMR = args.search_VF
path_to_AMR_database = args.VF_database

output_dir = args.output_dir

nthreads = args.threads

kraken2_col_names = ["percent_reads", "reads_rooted", "reads_taxon", "rank", "tax_id", "tax_name"]

bacterial_phyla = ['Bacillota','Actinomycetota','Bacteroidota','Pseudomonadota','Fusobacteriota','Campylobacterota']
dangerous_species = ['Bordetella pertussis', 'Burkholderia pseudomallei', 'Chlamydophila pneumoniae', 'Corynebacterium diphtheriae', 'Haemophilus influenzae', 'Mycoplasma pneumoniae', 'Streptococcus pneumoniae', \
                    'Legionella pneumophila', 'Staphylococcus aureus', 'Klebsiella pneumoniae', 'Moraxella catarrhalis', 'Pseudomonas aeruginosa', 'Acinetobacter baumannii', 'Stenotrophomonas maltophilia']



def format_kraken2_report(file_path):
    os.system('sed -E "s/^\s+//g" ' + file_path + '>' + file_path.split('.')[0] + '.tsv')
    os.system('sed -E -i "s/\t\s+/\t/g" ' + file_path.split('.')[0] + '.tsv')


def get_bacterial_info(dir_path, sample):

    df = pd.read_csv(os.path.join(dir_path, sample, sample + '.tsv'), sep="\t", names=kraken2_col_names)
    df = df[~df['rank'].isin(['S1', 'S2'])]
    df.loc[df['rank']=='S', 'reads_taxon'] = df.loc[df['rank']=='S', 'reads_rooted']
    df.to_csv(os.path.join(dir_path, sample, sample + '_species_filtered.tsv'), sep = '\t', header=False, index=False)

    reads_bacteria = df[df['tax_name']=='Bacteria']['reads_rooted'].tolist()[0] if not df[df['tax_name']=='Bacteria'].empty else 0
    perc_bacteria = df[df['tax_name']=='Bacteria']['percent_reads'].tolist()[0] if not df[df['tax_name']=='Bacteria'].empty else 0

    perc_by_phylum = [float(df[df['tax_name']==p]['reads_rooted'].tolist()[0])/reads_bacteria if not df[df['tax_name']==p].empty else 0 for p in bacterial_phyla]
    
    pathogens = ','.join(df[df['tax_name'].isin(dangerous_species) & (df['reads_rooted'] > 5000)]['tax_name'].tolist())

    df = df[~df['rank'].isin(['G1', 'G2', 'S'])]
    df.loc[df['rank']=='G', 'reads_taxon'] = df.loc[df['rank']=='G', 'reads_rooted']
    df.to_csv(os.path.join(dir_path, sample, sample + '_genus_filtered.tsv'), sep = '\t', header=False, index=False)

    return {"sample": sample, "reads_bacteria": reads_bacteria, "perc_bacteria": perc_bacteria, 'Bacillota (Firmicutes)': perc_by_phylum[0], 'Actinomycetota': perc_by_phylum[1],\
    'Bacteroidota': perc_by_phylum[2], 'Pseudomonadota (Proteobacteria)': perc_by_phylum[3], 'Campylobacterota': perc_by_phylum[4], 'Fusobacteriota': perc_by_phylum[5], \
    "pathogenic_species_threshold_reads": 5000, "pathogenic_species": pathogens}

def get_AMR_genes(dir_path, db_path, sample):
    os.system('ariba run ' +  db_path + ' ' + os.path.join(dir_path, sample, sample + '_merged_R1.fastq.gz') + ' ' + os.path.join(dir_path, sample, sample + '_merged_R2.fastq.gz') + ' ' +\
     os.path.join(dir_path, sample, 'ariba') + ' --threads ' + str(nthreads) + ' >/dev/null')


if __name__ == '__main__':

    df = pd.DataFrame(columns=['sample', 'reads_bacteria', 'perc_bacteria', 'Bacillota (Firmicutes)','Actinomycetota','Bacteroidota','Pseudomonadota (Proteobacteria)','Fusobacteriota','Campylobacterota', \
        'pathogenic_species_threshold_reads', 'pathogenic_species'])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    i = 0
    no_samples = len(os.listdir(path_to_kraken2_results))
    
    print('Processing samples:')
    for sample in os.listdir(path_to_kraken2_results):
        i += 1
        print(sample + '\t' + str(i) + '/' + str(no_samples))
        format_kraken2_report(os.path.join(path_to_kraken2_results, sample, sample + '.report'))
        if search_AMR:
            get_AMR_genes(path_to_raw_reads, path_to_AMR_database, sample)
        df = df._append(get_bacterial_info(path_to_kraken2_results, sample), ignore_index=True)

    meta = pd.read_csv(metadata_file, header = 0)
    df = meta.merge(df, on = "sample")
    
    if search_AMR:
        os.system('ariba summary ' + output_dir + '/ariba.summary ' + path_to_raw_reads + '/*/ariba/report.tsv --cluster_cols match,assembled --no_tree')

        amr = pd.read_csv(output_dir + '/ariba.summary.csv', header = 0)

        amr['name'] = [s.split('/')[1] for s in amr['name']]

        amr.columns = ['sample'] + [s.split('.')[0] for s in amr.columns[1:]]

        amr['virulence_factors'] = [','.join(amr.columns[1:][row[1:] != 'no']) for index, row in amr.iterrows()]
        amr['virulence_factor_count'] = [sum(row[1:-1] != 'no') for index, row in amr.iterrows()]

        amr = amr[['sample', 'virulence_factors', 'virulence_factor_count']]
        df = df.merge(amr, on='sample')
        df['bacterial_infection_criteria'] = 'min. 1 pathogenic species & min. 3 VFs'
        df['likely_bacterial_infection'] = ~(df['pathogenic_species'] == "") & (df['virulence_factor_count']>=3)


    df.to_csv(output_dir+'/plate_report.csv', index = False)
    