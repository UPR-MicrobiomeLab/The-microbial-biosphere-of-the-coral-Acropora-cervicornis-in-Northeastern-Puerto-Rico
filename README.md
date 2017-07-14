# The-microbial-biosphere-of-the-coral-Acropora-cervicornis-in-Northeastern-Puerto-Rico
Tutorial for analyzing the prokaryotic biomes of corals sampled at two depths in Puerto Rico.


* **TUTORIAL**

# A total of 1,594,650 reads with Phred scores>20 (in file seq_fna.fasta) underwent chimera identification using usearch61, which will provide a list of chimeric sequences. 
You can download these sequences from NCBI BioProject ID PRJNA379103. 

identify_chimeric_seqs.py -m usearch61 -i seq_fna.fasta --suppress_usearch61_ref -o chimera_checking

# Chimeras were then removed using the “ -n” flag

filter_fasta.py -f seq_fna.fasta -o seqs_chimeras_filtered_corales.fna -s chimera_checking/chimeras.txt -n

# Non-chimeric reads (n=1,538,471 in file seqs_chimeras_filtered_corales.fna from previous script) underwent OTU picking with usearch61. 

pick_otus.py -m usearch61 -i seqs_chimeras_filtered_corales.fna -o usearch61_picked_otus

# Representative set of sequences per OTU was selected using script pick_rep_set.py. The resulting rep_set.fna file is available from this github project.

pick_rep_set.py -i usearch61_picked_otus/seqs_chimeras_filtered_otus.txt -f split_library/seq_fna.fasta -o rep_set.fna

# Greengenes Taxonomy was assigned using the August 2013 Greengenes database with the putput being the full lineage information for each representative sequence 

assign_taxonomy.py -i rep_set.fna -r Gg_13_8_99.taxonomy/gg_13_8_99.fasta -t Gg_13_8_99.taxonomy/gg_13_8_99.gg.tax 

# Sequences greater then 200bp were aligned with 97% similarity.

align_seqs.py -i rep_set.fna -e 200 -p 97.0

#The alignment was filtered to remove positions which are gaps in every sequence, to differentiate between non-conserved positions, which are uninformative for tree building, and conserved positions which are informative for tree building.

filter_alignment.py -i pynast_aligned/rep_set_aligned.fasta -o filtered_alignment

# A phylogenetic tree was built. The resulting rep_set.tree file is available from this tutorial. 

make_phylogeny.py -i filtered_alignment/rep_set_aligned_pfiltered.fasta -o rep_set.tre

# We then proceeded to build an OTU table with the passing sequences from the previous aligned and taxonomy assigned steps. The otu_table-inicial.biom is also provided in this tutorial.

make_otu_table.py -i usearch61_picked_otus/seqs_chimeras_filtered_otus.txt -t uclust_assigned_taxonomy/rep_set_tax_assignments.txt -e pynast_aligned/rep_set_failures.fasta -o otu_table-inicial.biom

#A total of 739,664 good-quality sequences were reflected in the OTU table, however these included Chloroplasts and singletons which may reflect sequencing artifacts, therefore we proceeded to eliminate OTUs manually after converting the biom to a text file. 

biom convert –i otu_table-inicial.biom –o otu_table-inicial.txt --to-tsv --header-key taxonomy

#We opened the txt file in excel and removed mannually: 1) the only OTU corresponding to an Archaea (present in only 1 sample), 2) kept OTUs with sequences present in at least three of the 6 samples. We additionally proceeded with further manual filtering removing removed OTUs with less than three sequences per sample. 
Such conservative analysis was intended to reduce the effect of possible sequencing artifacts.
It yielded a conserved and simplified OTU table (that mainly removed low count Ricketsialles OTUs), resulting in 803 OTUs with 173,137 sequences that were used in further analyses that was saved as "coral.txt". This file is also available in this tutorial as is its biom version "coral.biom" used in further analyses.

#The file was saved as coral.txt and converted to coral.biom using the script:

biom convert –i coral.txt –o coral.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#Interestingly we found populations of Ricketsialles-dominant taxa and other non-Ricketsialles taxa. We proceeded with the data analyses jointly as well as separating the rare from the dominant taxa. We have provided an excel file Coral-PR-Otutables.xlsx that summarizes the OTU table used in the analyses as well as separated the and the rare and ricketsialles-dominant taxa (in separate sheets). 

#From now on we will detail the data analyses used to generate the figures in this manuscript:
 
#Number of sequences and OTUs used in the analyses as reflected in Table 1.
#To evaluate the statistics of the OTU table regarding the number of sequences per sample, we used:

biom summarize-table -i  coral.biom

# To evaluate the OTU statistics we used the same script with the qualitative flag:

biom summarize-table -i coral.biom --qualitative



