#

configfile: "config.json"

#######################################################################################
#   Prepare human references
#######################################################################################

rule prepare_reference_human_caNFATC1:
    output:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='caNFATC1'
        )
    run:
        commands = [
            'mkdir -p {out}/{viral}/rsem_reference/ref'.format(
                out=config['output'],
                viral='caNFATC1'
            ),
            'rsem-prepare-reference --bowtie {ref},{tran}/{viral}/{viral}.fa {out}/{viral}/rsem_reference/ref'.format(
                ref=config['human_transcriptome'],
                tran=config['transcript_files'],
                viral='caNFATC1',
                out=config['output']
            )
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)

rule prepare_reference_human_caNFATC2v1:
    output:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='caNFATC2_v1'
        )
    run:
        commands = [
            'mkdir -p {out}/{viral}/rsem_reference/ref'.format(
                out=config['output'],
                viral='caNFATC2_v1'
            ),
            'rsem-prepare-reference --bowtie {ref},{tran}/{viral}/{viral}.fa {out}/{viral}/rsem_reference/ref'.format(
                ref=config['human_transcriptome'],
                tran=config['transcript_files'],
                viral='caNFATC2_v1',
                out=config['output']
            )
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)

rule prepare_reference_human_no_viral:
    output:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='no_viral'
        )
    run:
        commands = [
            'mkdir -p {out}/{viral}/rsem_reference/ref'.format(
                out=config['output'],
                viral='no_viral'
            ),
            'rsem-prepare-reference --bowtie {ref} {out}/{viral}/rsem_reference/ref'.format(
                ref=config['human_transcriptome'],
                tran=config['transcript_files'],
                viral='no_viral',
                out=config['output']
            )
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)


###################################################################################
#   Quantify human samples
###################################################################################
rule quantify_human_caNFATC1:
    input:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='caNFATC1'
        )
    output:
        expand(
            '{out}/{viral}/quantification/{{prefix}}.isoforms.results'.format(
                out=config['output'],
                viral='caNFATC1'
            ),
            prefix=config['human_caNFATC1_prefixes']
        )
    run:
        commands=['mkdir -p {}/caNFATC1/quantification'.format(config['output'])]
        commands+=[
            'nohup rsem-calculate-expression {data}/{prefix}.fastq  {out}/caNFATC1/rsem_reference/ref {out}/{viral}/quantification/{prefix} &'.format( 
                viral='caNFATC1',
                data=config['human_raw_data'],
                out=config['output'],
                prefix=prefix
            )
            for prefix in config['human_caNFATC1_prefixes']
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)

rule quantify_human_caNFATC2_v1:
    input:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='caNFATC2_v1'
        )
    output:
        expand(
            '{out}/{viral}/quantification/{{prefix}}.isoforms.results'.format(
                out=config['output'],
                viral='caNFATC2_v1'
            ),
            prefix=config['human_caNFATC2_prefixes']
        )
    run:
        commands=['mkdir -p {}/caNFATC2_v1/quantification'.format(config['output'])]
        commands+=[
            'nohup rsem-calculate-expression {data}/{prefix}.fastq  {out}/caNFATC2_v1/rsem_reference/ref {out}/{viral}/quantification/{prefix} &'.format(
                viral='caNFATC2_v1',
                data=config['human_raw_data'],
                out=config['output'],
                prefix=prefix
            )
            for prefix in config['human_caNFATC2_prefixes']
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)


rule quantify_human_no_viral:
    input:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='human_NFATC1',
            viral='no_viral'
        )
    output:
        expand(
            '{out}/{viral}/quantification/{{prefix}}.isoforms.results'.format(
                out=config['output'],
                viral='no_viral'
            ),
            prefix= config['human_caNFATC1_prefixes'] + config['human_caNFATC2_prefixes'] + config['human_GFP_prefixes']
        )
    run:
        commands=['mkdir -p {}/no_viral/quantification'.format(config['output'])]
        commands+=[
            'nohup rsem-calculate-expression {data}/{prefix}.fastq  {out}/{viral}/rsem_reference/ref {out}/{viral}/quantification/{prefix} &'.format(
                viral='no_viral',
                data=config['human_raw_data'],
                out=config['output'],
                prefix=prefix
            )
            for prefix in config['human_caNFATC1_prefixes'] + config['human_caNFATC2_prefixes'] + config['human_GFP_prefixes']
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)


#######################################################################################
#   Prepare mouse references
#######################################################################################

rule prepare_reference_mouse_caNFATC1:
    output:
        '{out}/{viral}/rsem_reference/ref.transcripts.fa'.format(
            out=config['output'],
            endo='mouse_NFATC1',
            viral='caNFATC1'
        )
    run:
        commands = [
            'mkdir -p {out}/{viral}/rsem_reference/ref'.format(
                out=config['output'],
                viral='caNFATC1'
            ),
            'rsem-prepare-reference --bowtie {ref},{tran}/{viral}/{viral}.fa {out}/{viral}/rsem_reference/ref'.format(
                ref=config['mouse_transcriptome'],
                tran=config['transcript_files'],
                viral='caNFATC1',
                out=config['output']
            )
        ]
        for c in commands:
            shell('echo "{}"'.format(c))
            shell(c)


