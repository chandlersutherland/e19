# e19
library prep for somatic mutation sequencing 

`re_simulation.py`: simulates digestion of genome, outputing a fasta file of fragments of the desired size, and a corresponding bed file. Also an option to simulate 150PE reads on these fragments, outputing the actual sequenced reads is available. 

`re_mapping.sh`: map simulated reads to the genome, and use samtools coverate to generate a .tsv file for the whole genome, the gene space, and NLRs specifically. 

`re_processing.py`: combine results of `re_mapping.sh` into a single file 

`re_plotting.Rmd`: generate some figures to look at the simulation results 

`abascal_enzymes.txt`: enzymes tested here and in Abascal 2021. 
