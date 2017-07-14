# The-microbial-biosphere-of-the-coral-Acropora-cervicornis-in-Northeastern-Puerto-Rico
Tutorial for analyzing the prokaryotic biomes of corals sampled at two depths in Puerto Rico.


* **TUTORIAL**

Data analyses includes standard python QIIME scripts as well as R packages such as VEGAN,ggplot2,reshape2, devtools, RColorBrewer. This tutorial contains the logical workflow, scripts and necessary intermediary files for the production of the data discussed in the manuscript titled: "The microbial biosphere of the coral _Acropora cervicornis_ in Northeastern Puerto Rico" by Godoy-Vitorino, F., Ruiz-Diaz,C.P., Rivera-Seda, A., Ramírez-Lugo, J.S. and Toledo-Hernández, C.

#A total of 1,594,650 reads with Phred scores>20 (in file seq_fna.fasta) underwent chimera identification using usearch61, which will provide a list of chimeric sequences. 
# Sequences can be downloaded NCBI SRA Archive BioProject ID PRJNA379103. 

# Identify Chimeras from your raw sequences with usearch61 in QIIME

identify_chimeric_seqs.py -m usearch61 -i seq_fna.fasta --suppress_usearch61_ref -o chimera_checking

# Remove Chimeras using the “ -n” flag in QIIME

filter_fasta.py -f seq_fna.fasta -o seqs_chimeras_filtered_corales.fna -s chimera_checking/chimeras.txt -n

# Pick OTUs in QIIME
#Non-chimeric reads (n=1,538,471 in file seqs_chimeras_filtered_corales.fna from previous script) underwent OTU picking with usearch61. 

pick_otus.py -m usearch61 -i seqs_chimeras_filtered_corales.fna -o usearch61_picked_otus

# Select a set of representative sequences per OTU
#Representative set of sequences per OTU was selected using script pick_rep_set.py. The resulting rep_set.fna file is available from this github project.

pick_rep_set.py -i usearch61_picked_otus/seqs_chimeras_filtered_otus.txt -f split_library/seq_fna.fasta -o rep_set.fna

# Assign taxonomy to your representative OTUs
#Greengenes Taxonomy was assigned using the August 2013 Greengenes database with the putput being the full lineage information for each representative sequence 

assign_taxonomy.py -i rep_set.fna -r Gg_13_8_99.taxonomy/gg_13_8_99.fasta -t Gg_13_8_99.taxonomy/gg_13_8_99.gg.tax 

# Align sequences
#Sequences greater then 200bp were aligned with 97% similarity.

align_seqs.py -i rep_set.fna -e 200 -p 97.0

# Filter the alignment to remove positions which are gaps in every sequence, to differentiate between non-conserved positions, which are uninformative for tree building, and conserved positions which are informative for tree building.

filter_alignment.py -i pynast_aligned/rep_set_aligned.fasta -o filtered_alignment

# Build a phylogenetic tree.
#The resulting rep_set.tree file is available from this tutorial. 

make_phylogeny.py -i filtered_alignment/rep_set_aligned_pfiltered.fasta -o rep_set.tre

# Build an OTU table with the passing sequences from the previous aligned and taxonomy assigned steps. 
#The otu_table-inicial.biom is also provided in this tutorial.

make_otu_table.py -i usearch61_picked_otus/seqs_chimeras_filtered_otus.txt -t uclust_assigned_taxonomy/rep_set_tax_assignments.txt -e pynast_aligned/rep_set_failures.fasta -o otu_table-inicial.biom

#A total of 739,664 good-quality sequences were reflected in the OTU table, however these included Chloroplasts and singletons which may reflect sequencing artifacts, therefore we proceeded to eliminate OTUs manually after converting the biom to a text file. 

# Convert your initial OTU table to a text format and do manual curation.

biom convert –i otu_table-inicial.biom –o otu_table-inicial.txt --to-tsv --header-key taxonomy

#We opened the txt file in excel and removed mannually: 1) the only OTU corresponding to an Archaea (present in only 1 sample), 2) kept OTUs with sequences present in at least three of the 6 samples. We additionally proceeded with further manual filtering removing removed OTUs with less than three sequences per sample. Such conservative analysis was intended to reduce the effect of possible sequencing artifacts.

# The resuting simplified OTU table resulted in 803 OTUs with 173,137 sequences that were used in further analyses that was saved as "coral.txt". 
#This file is also available in this tutorial as is its biom version "coral.biom" used in further analyses.
#The file was saved as coral.txt and converted to coral.biom using the script:

biom convert –i coral.txt –o coral.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#Interestingly we found populations of Ricketsialles-dominant taxa and other non-Ricketsialles taxa. We proceeded with the data analyses jointly as well as separating the rare from the dominant taxa. 
We have provided an excel file Coral-PR-Otutables.xlsx that summarizes the OTU table used in the analyses as well as separated the and the rare and ricketsialles-dominant taxa (in separate sheets). 

# From now on we will detail the data analyses used to generate the figures in this manuscript
 
#Number of sequences and OTUs used in the analyses as reflected in Table 1.
# To evaluate the statistics of the OTU table regarding the number of sequences per sample, we used:

biom summarize-table -i  coral.biom

# To evaluate the OTU statistics we used the same script with the qualitative flag:

biom summarize-table -i coral.biom --qualitative

# For Figure 2 
#Figure 2 beta and alpha div.pdf) available in this repository.
#Files used: 1) coral.biom; 2) rep_set.tre; 3) Coral_Mapping.txt; 4) alpha_params.txt

# For Fig 2 Panel C: 
#Alpha rarefaction curves in QIIME followed by statistical testing

alpha_rarefaction.py -i coral.biom -p alpha_params.txt -t rep_set.tre -m Coral_Mapping.txt -o rare/ -e 28000

compare_alpha_diversity.py -i rare/alpha_div_collated/PD_whole_tree.txt -m Coral_Mapping.txt -c Depth_m -d 28000 -o statpd

#The Depth_m_stats.txt file is available in this repository.

# Figure 2B NMDS:
#files needed are: 1) Distance matrix (38distmatrix.csv). Built through the L6 genus-level taxa table (of the 78 non-ricketsialles OTUs) converted to text with taxonomy in the first column 2) Metadata file: metacoral.txt
#Run nmds-for-fig2B.R (also available here) with the above files.

# For figure 2C 3D-Biplots:
#These are made with emperor in QIIME files 1) norick.biom = 78 non-ricketsialles OTUs; rep_set.tre;  Coral_Mapping.txt

beta_diversity_through_plots.py -i norick.biom -t rep_set.tre -m Coral_Mapping.txt -o beta

make_emperor.py -i beta/weighted_unifrac_pc.txt -m Coral_Mapping.txt -o biplotoptions6 -t taxarare/L6.txt -n 6 --add_vectors depth

# Figure 3 Taxa Summaries of filtered Rickettsiales and Non-filtered ricketsialles taxa: 
#File rarefied_obs_table.biom; norick_obs_table.biom available here will be used to generate taxa summary plots in QIIME that will have L2 and L6 table used in R to generate the figures.

single_rarefaction.py -i coral.biom -d 28000 -o rarefied_obs_table.biom

summarize_taxa_through_plots.py -i rarefied_obs_table.biom -m Coral_Mapping.txt -o rarifiedtaxatsum 

filter_taxa_from_otu_table.py -i rarefied_obs_table.biom -n o__Rickettsiales -o norick_obs_table.biom

summarize_taxa_through_plots.py -i norick_obs_table.biom -m Coral_Mapping.txt -o taxa_summary_norick

#To visualize taxa summary bar plots we used GGPLOT2 R package. The data was prepared for plotting by 'melting' using the RESHAPE2 R package. We merged the L2 and filtered Rickettsiales L2 output from QIIME for the phyla plots and we merged the L6 and filtered Rickettsiales L6 for the genus plots. 
#This melting was done by columns Taxonomy, the variable for color, and by FACET (to establish two panels), the variable for splitting into facets, or panels, to highlight non-Rickettsiales organisms. Color palette was done manually with the assistance of colorpicker.com. 
#The files needed to generate Fig 3 are called “taxonomy_table_L6_R” and “taxonomy_table_L2_R” also available here as well as the scripts used to generate the files, taxa_summary_phyla_script.R and taxa_summary_genus_script.R.

# Figure 4 - Taxa summary of the OTUs classified as Rickettsiales: 
#Files needed and provided here: 1)ricketsialles.txt; 2) ricketsialles.biom; 3) rickrarefied.biom; 4) gtest_depth-ricketsialles-44; 5) sigrick44.biom and Figure generated “Figure 4-ricketsialles.pdf”

biom convert -i ricketsialles.txt -o ricketsialles.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy 

single_rarefaction.py -i ricketsialles.biom -d 23472 -o rickrarefied.biom

summarize_taxa_through_plots.py -i rickrarefied.biom -m Coral_Mapping.txt -o taxarick/

alpha_rarefaction.py -i rickrarefied.biom -p alpha_params.txt -t rep_set.tre -m Coral_Mapping.txt -o rarerick/ -e 23472

group_significance.py -i rickrarefied.biom -m Coral_Mapping.txt -c Depth_m -s g_test -o gtest_depth-ricketsialles.txt

filter_otus_from_otu_table.py -i rickrarefied.biom -e gtest_depth-ricketsialles-44.txt -o sigrick44.biom --negate_ids_to_exclude

normalize_table.py -i sigrick44.biom -a DESeq2 -o DESeq2_g_test44.biom

biom convert -i DESeq2_g_test44.biom -o DESeq2_g_test44.txt --to-tsv --header-key taxonomy

#File DESeq2_g_test44_tax.txt will be used to prepare the CSV file for Heatmap.3 plotting in R. 1st run file heatmap_function.R and after, apply heatmap-script-fig4.R both available in this repository.


# For Figure 5 Core OTUs:
#you will need to use file norick_obs_table.biom.

compute_core_microbiome.py -i norick_obs_table.biom -o core-deep/ --mapping_fp Coral_Mapping.txt --valid_states "depth:deep"

compute_core_microbiome.py -i norick_obs_table.biom -o core-shallow/ --mapping_fp Coral_Mapping.txt --valid_states "depth:shallow"
 
#Manually merge OTU IDs from core_otus_100s.txt in both depths. Save the list of core OTU IDs-coretaxa “list-coretaxa.TXT” and create a new OUT table with only these specific core OTUs

biom convert -i rarefied_obs_table.txt -o rarefied_obstable.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

filter_otus_from_otu_table.py -i norick_rarefied_obs_table.biom -e list-coretaxa.txt -o corecoral.biom --negate_ids_to_exclude

summarize_taxa_through_plots.py -i corecoral.biom -m Coral_Mapping.txt -o taxacorecoral/

# For Figure 6 - Panel A:
#Heatmap showing the 38 significantly different taxa between shallow and deep samples Panel B shows selected boxplots highlighting different taxa that changes significantly according to the coral depth.

group_significance.py -i norick_rarefied_obs_table.biom -m Coral_Mapping.txt -c Depth_m -s g_test -o gtest_depth.txt

#Remember the selected p-values were FDR corrected <0.01

filter_otus_from_otu_table.py -i norick_rarefied_obs_table.biom -e gtest_depth.txt -o signorickrick38.biom --negate_ids_to_exclude

biom convert -i signorickrick38.biom -o signorickrick38.txt --to-tsv --header-key taxonomy

normalize_table.py -i signorickrick38.biom -a DESeq2 -o DESeq2_g_test38.biom

biom convert -i DESeq2_g_test38.biom -o DESeq2_g_test38.txt --to-tsv --header-key taxonomy

#note: Complete the DESeq2_g_test38.txt taxonomy column with the signorickrick38.txt info. Manually Prepare the CSV file for R retaining only genus, resulting file: gtest38.csv and use Coral_Mapping.csv

#1st run file heatmap_function.R  and follow with the implementation of script heatmap-script-fig6.R. All files needed are available in this repository.

# For Figure 6 - Panel B (boxplots): 
#an example is give below. Use the gtest38.csv and select lines for the taxa you’re interested in (any OTU from the heatmap in taxa can work), just follow example given here is _Serratia marcescens_, save the row as a text file serratia-marcescens.txt FDR p-values are written according to their values in the gtest_depth.txt file provided before. Use library vegan and run script provided in Boxplots-sig-dif-taxa.R.

# Congratulations!!



# LIST OF ALL THE 45 FILES AVAILABLE IN THIS REPOSITORY (SCRIPTS, FIGURES AND INTERMEDIARY FILES):
	README.md
	Figure 2 beta and alpha div.pdf
	Figure 4-ricketsialles.pdf
	Figure 5 core 100.pdf
	Figure3-taxaplots.pdf
	Figure6-small.tif
	38distmatrix.csv
	Boxplots-sig-dif-taxa.R
	Coral-PR-Otutables.xlsx
	Coral_Mapping.csv
	Coral_Mapping.txt
	DESeq2_g_test38.biom
	DESeq2_g_test38.txt
	DESeq2_g_test44.txt
	DESeq2_g_test44_tax.csv
	Depth_m_stats.txt
	alpha_params.txt
	coral.biom
	coral.txt
	gtest38.csv
	gtest_depth-ricketsialles-44.txt
	gtest_depth.txt
	heatmap-script-fig4.R
	heatmap-script-fig6.R
	heatmap_function.R
	list-coretaxa.txt
	metacoral.txt
	nmds-for-fig2B.R
	norick.biom
	norick_obs_table.biom
	norick_rarefied_obs_table.biom
	otu_table-inicial.biom
	rarefied_obs_table.biom
	rep_set.fna
	rep_set.tre
	ricketsialles.biom
	ricketsialles.txt
	rickrarefied.biom
	signorickrick38.biom
	signorickrick38.txt
	sigrick44.biom
	taxa_summary_genus_script.R
	taxa_summary_phyla_script.R
	taxonomy_L2_R.txt
	taxonomy_L6_R.txt

