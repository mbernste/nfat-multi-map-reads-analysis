This code runs RSEM on all of the human and mouse raw RNA-seq islet data using references that include the inserted viral transcript. To re-run this analysis, install the Python packages in `requirements.txt`:  

`pip install -r requirements.txt`  

This analysis requires RSEM v1.2.21 and Bowtie v1.2.3 are installed and accessible via the `PATH` variable. 

This pipeline uses [Snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) to coordinate all of the building of the references and running RSEM. To run the full pipeline from beginning to end, we will need to first configure the pipeline to point the reference data, raw data, and to where we want to write the output. This is all specified in the `config.json` file. The following fields will need to be adjusted:  
- human_transcriptome: The path to the human transcriptome FASTA file  
- mouse_transcriptome: The path to the mouse transcriptome FASTA file  
- human_raw_data: The path to the human raw reads FASTQ files  
- mouse_raw_data: The path to the mouse raw reads FASTQ files  
- output: The path to where all output should be written

To run the pipeline, run the following command from the project directory:  

`snakemake`  

Note, this will take some time to complete as it builds the reference, runs the alignment, and quantifies expression for 39 samples.  

The DAG depicting the workflow:   

![screenshot](https://github.com/mbernste/nfat-multi-map-reads-analysis/blob/master/dag.png)

