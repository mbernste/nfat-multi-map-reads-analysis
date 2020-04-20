This code runs RSEM on all of the human and mouse raw RNA-seq islet data using references that include the inserted viral transcript. To re-run this analysis, install the Python packages in `requirements.txt`:  
`pip install -r requirements.txt`  
This analysis requires RSEM v1.2.21 and Bowtie v1.2.3 are installed and accessible via the `PATH` variable. 

This pipeline uses [Snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) to coordinate all of the building of the references and running RSEM. To run the full pipeline from beginning to end, run   
`snakemake`  
in the project's root directory. Note, this will take a long time to complete as it builds the reference, runs the alignment, and quantification for 39 samples.   
![screenshot](https://github.com/mbernste/nfat-multi-map-reads-analysis/blob/master/dag.png)

