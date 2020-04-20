import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import pysam
from collections import defaultdict, Counter
import upsetplot
from upsetplot import from_memberships
import sys
import json
from os.path import join
import pandas as pd

def main():
    with open('config.json', 'r') as f:
        config = json.load(f)

    human_endo_genes = [
        'human_NFATC1',
        'human_NFATC2'
    ]
    mouse_endo_genes = [
        'mouse_NFATC1',
        'mouse_NFATC2'
    ]
    human_canfatc1_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/caNFATC1/quantification' #args[1]
    human_canfatc2_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/caNFATC2_v1/quantification'
    human_no_viral_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/no_viral/quantification' #args[2]
    mouse_canfatc1_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/mouse_caNFATC1/quantification' #args[1]
    mouse_canfatc2_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/mouse_caNFATC2_v1/quantification'
    mouse_no_viral_dir = '/scratch/mnbernstein/attie_collab/hg19_alignment/mouse_no_viral/quantification' #args[2]
    out_dir = '.'

    human_canfatc1_df = _gather_results(
        human_endo_genes, 
        'caNFATC1', 
        human_canfatc1_dir, 
        human_no_viral_dir,
        config['human_caNFATC1_prefixes'],
        config['transcript_metadata']
    )
    human_canfatc2_df = _gather_results(
        human_endo_genes, 
        'caNFATC2_v1', 
        human_canfatc2_dir, 
        human_no_viral_dir,
        config['human_caNFATC2_prefixes'],
        config['transcript_metadata']
    )
    mouse_canfatc1_df = _gather_results(
        mouse_endo_genes,
        'caNFATC1',
        mouse_canfatc1_dir,
        mouse_no_viral_dir,
        config['mouse_caNFATC1_prefixes'],
        config['transcript_metadata']
    )

    fig, axarr = plt.subplots(2,2,figsize=(9,8))
    sns.set(style="whitegrid", palette="pastel")
    sns.set_style("whitegrid")
    sns.boxplot(
        x="Gene", 
        y="TPM", 
        hue="Reference", 
        palette=["m", "g"],
        ax=axarr[0][0],
        data=human_canfatc1_df
    )
    axarr[0][0].get_legend().remove()
    axarr[0][0].set_yscale('log')
    axarr[0][0].set_title('Human caNFATC1-treated')

    sns.boxplot(
        x="Gene",
        y="TPM",
        hue="Reference",
        palette=["m", "g"],
        ax=axarr[0][1],
        data=human_canfatc2_df
    )
    axarr[0][1].set_yscale('log')
    axarr[0][1].set_title('Human caNFATC2-treated')
    axarr[0][1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    sns.boxplot(
        x="Gene",
        y="TPM",
        hue="Reference",
        palette=["m", "g"],
        ax=axarr[1][0],
        data=mouse_canfatc1_df
    )
    axarr[1][0].set_yscale('log')
    axarr[1][0].set_title('Mouse caNFATC1-treated')
    axarr[1][0].get_legend().remove()

    #axarr[0][0].grid(b=True, which='minor', axis='y')

    plt.tight_layout()
    fig.savefig(
        join(out_dir, 'include_viral_vs_no_include_viral.png'),
        format='png',
        dpi=150
    )


def _gather_results(endo_genes, viral_gene, viral_dir, no_viral_dir, prefixes, gene_to_transcripts):
    da = []
    for prefix in prefixes:
        print('Reading expression data for sample ', prefix)
        for gene in endo_genes:
            viral_expr = _retrieve_expression(
                join(viral_dir, '{}.isoforms.results'.format(prefix)),
                gene_to_transcripts[gene]
            )
            no_viral_expr = _retrieve_expression(
                join(no_viral_dir, '{}.isoforms.results'.format(prefix)),
                gene_to_transcripts[gene]
            )
            da.append((
                prefix, 'Include viral', gene.split('_')[1], viral_expr
            ))
            da.append((
                prefix, 'Exclude viral', gene.split('_')[1], no_viral_expr
            ))
    df = pd.DataFrame(
        data=da,
        columns=[
            'Sample',
            'Reference',
            'Gene',
            'TPM'
        ]
    )
    print(df)
    return df


def _retrieve_expression(fname, isoforms):
    df = pd.read_csv(fname, sep='\t', index_col=0)
    return sum(df.loc[isoforms]['TPM'])
   

if __name__ == "__main__":
    main() 
