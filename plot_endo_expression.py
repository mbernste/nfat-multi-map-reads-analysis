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

    ENDO_GENES = [
        'human_NFATC1',
        'human_NFATC2'
    ]
    viral_quants = '/scratch/mnbernstein/attie_collab/hg19_alignment/caNFATC1/quantification' #args[1]
    no_viral_quants = '/scratch/mnbernstein/attie_collab/hg19_alignment/no_viral/quantification' #args[2]
    viral_gene = 'caNFATC1'
    out_dir = '.'

    da = []
    for prefix in config['human_caNFATC1_prefixes']:
        print('Reading expression data for sample ', prefix)
        for gene in ENDO_GENES:
            viral_expr = _retrieve_expression(
                join(viral_quants, '{}.isoforms.results'.format(prefix)),
                config['transcript_metadata'][gene]
            )
            no_viral_expr = _retrieve_expression(
                join(no_viral_quants, '{}.isoforms.results'.format(prefix)),
                config['transcript_metadata'][gene]
            )
            da.append((
                prefix, 'Include {}'.format(viral_gene), gene, viral_expr  
            ))
            da.append((
                prefix, 'Exclude {}'.format(viral_gene), gene, no_viral_expr  
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

    fig, ax = plt.subplots(1,1,figsize=(6,4))
    sns.set(style="ticks", palette="pastel")
    sns.boxplot(
        x="Gene", 
        y="TPM", 
        hue="Reference", 
        palette=["m", "g"],
        data=df
    )
    fig.get_axes()[0].set_yscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title('Expression for {}-treated samples'.format(viral_gene))
    plt.tight_layout()
    fig.savefig(
        join(out_dir, 'caNFATC1_include_viral_vs_no_include_viral.png'),
        format='png',
        dpi=150
    )




def _retrieve_expression(fname, isoforms):
    df = pd.read_csv(fname, sep='\t', index_col=0)
    return sum(df.loc[isoforms]['TPM'])
   

if __name__ == "__main__":
    main() 
