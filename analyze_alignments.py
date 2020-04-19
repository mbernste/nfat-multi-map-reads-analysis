import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import pysam
from collections import defaultdict, Counter
import upsetplot
from upsetplot import from_memberships
import sys
import json
from os.path import join

with open('config.json', 'r') as f:
    config = json.load(f)

endo_genes = sys.argv[1].split(',')
viral_gene = sys.argv[2]
bam_f = sys.argv[3]
out_pref = sys.argv[4]

endo_trans = []
for endo_gene in endo_genes:
    endo_trans += config['transcript_metadata'][endo_gene]
viral_trans = config['transcript_metadata'][viral_gene]

print('Reading alignment from ', bam_f) 
samfile = pysam.AlignmentFile(
    bam_f, "rb"
)
all_transcripts = endo_trans + viral_trans

read_to_transcripts = defaultdict(lambda: set())
for gene in endo_genes:
    for tran in config['transcript_metadata'][gene]:
        print('Processing reads aligned to {}...'.format(tran))
        reads = [r.to_dict()['name'] for r in samfile.fetch(tran)]
        for read in reads:
            read_to_transcripts[read].add(gene)
for tran in viral_trans:
    print('Processing reads aligned to {}...'.format(tran))
    reads = [r.to_dict()['name'] for r in samfile.fetch(tran)]
    for read in reads:
        read_to_transcripts[read].add(viral_gene)

memb_sets = [
    frozenset(v)
    for v in read_to_transcripts.values()
]
memb_set_count = Counter(memb_sets)
#print(memb_set_count)

#if len(memb_set_count.keys()) > 20:
#    print('Filtering reads for those mapping to small transcript intersections...')
#    filt = set([
#        frozenset(memb_set)
#        for memb_set, count in memb_set_count.items()
#        if count < 10
#    ])
#    plot_read_to_transcripts = {
#        read: transcripts
#        for read, transcripts in read_to_transcripts.items()
#        if frozenset(transcripts) not in filt 
#    }
#    print('done.')
#else:
plot_read_to_transcripts = read_to_transcripts    

memb = from_memberships(plot_read_to_transcripts.values())
upsetplot.plot(memb, subset_size='count', show_counts=True, sort_by='cardinality')
out_f = '{}.upset.png'.format(out_pref)
print('Plotting Upset to ', out_f)
plt.savefig(
    out_f, 
    format='png'
)

n_multimap_viral_endo = sum([
    1
    for transcripts in read_to_transcripts.values()
    if viral_gene in set(transcripts)
    and len(transcripts) > 1
])
n_mapped_endo = sum([
    1
    for transcripts in read_to_transcripts.values()
    if len(set(transcripts) & set(endo_trans)) > 1
])
n_mapped_viral = sum([
    1
    for transcripts in read_to_transcripts.values()
    if len(set(transcripts) & set(viral_trans)) > 1
])
out_f = '{}.report.txt'.format(out_pref)
print('Writing to ', out_f)
with open(out_f, 'w') as f:
    json.dump(
        {
            'multi_map_{}_{}'.format(endo_gene, viral_gene): n_multimap_viral_endo,
            'total_mapped': len(read_to_transcripts),
            'total_mapped_endogenous': n_mapped_endo,
            'total_mapped_viral': n_mapped_viral
        },
        f,
        indent=True
    )
