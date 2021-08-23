# Breaking the restriction barriers and applying CRISPRi as a gene silencing tool in Pseudoclostridium thermosuccinogenes

#### Joyshree Ganguly<sup>1</sup>, Maria Martin-Pascual<sup>2</sup>, Diego Montiel González<sup>3</sup>, Alkan Bulut<sup>3</sup>, Bram Vermeulen<sup>2</sup>, Ivo Tjalma<sup>2</sup>, Athina Vidaki<sup>3,c</sup>, Richard van Kranenburg<sup>1,2</sup>

* 1 Corbion, Arkelsedijk 46, 4206 AC Gorinchem, The Netherlands 
* 2 Laboratory of Microbiology, Wageningen University and Research, Stippeneng 4, 6708 WE Wageningen, The Netherlands
* 3 Department of Genetic Identification, Faculty of Medicine, Erasmus University Rotterdam, 3000 CB Rotterdam, The Netherlands
* 4 Fontys University of Applied Sciences, 5612 AR Eindhoven, The Netherlands
* contributed equally 


## Installation requirements 


> BLAST+ v2.6.0 <br>
> python>=3.6 & anaconda3  <br>
> R language >=3.6.1 "Action of the toes"
    

## Scripts and pipieline


1) Obtain list of motifs combinations based on the substitution DNA criteria    
    Usage: 
    <br>
    
    > python get_motifs_combinations.py
    
    Output: <br>
    > Files with a list of motif combinations 

    * [Included in motifs folder with extension .txt]

2) 2.1 Blast: Create a database index reference for e coli with following command:

    > makeblastdb -in ecoli/HST04_SEQt.fasta -out ecoli/HST04 -dbtype nucl 

   2.2 Performs a blastn search with the motifs combinations against the ecoli reference:
    
	> blastn -query motifs/CACNNNNNNNTNGC.txt -db ecoli/HST04 -out motifs/CACNNNNNNNTNGC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 4;
    
	> blastn -query motifs/DGAGNNNNATC.txt -db ecoli/HST04 -out motifs/DGAGNNNNATC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 4;
	blastn -query motifs/GCNANNNNNNNGTG.txt -db ecoli/HST04 -out motifs/GCNANNNNNNNGTG.out -task blastn-short -outfmt 7 -perc_identity 100 -num_threads 4;
    
	> blastn -query motifs/TCABNNNNNNTARG.txt -db ecoli/HST04 -out motifs/TCABNNNNNNTARG.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 4;

	> blastn -query motifs/CYTANNNNNNVTGA.txt -db ecoli/HST04 -out motifs/CYTANNNNNNVTGA.out -task blastn-short -outfmt 7 -perc_identity 100 -num_threads 4;

	> blastn -query motifs/GATNNNNCTC.txt -db ecoli/HST04 -out motifs/GATNNNNCTC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 4;

    Output: 
    >Table of hits for each motifs list (combinations) against the ecoli reference database

    * [included in motifs folder with extension .out]

3) Get A modification based on list of motifs combinations and wig files
    Usage:
    > python get_modifications.py motifs/ 0 4

    >motifs/: is the directory containing motifs combinations prodcued from get_motifs_combinations script <br>
    mismatches [0]: number of mismatches allowed at the end of the motif     <br>
    n_threads[4]: number of cpus to use during get modifications 
    * [This can be very exhaustive by looking at almost 1,000,000 (one million) positions of each wig to 200,000 motifs combinations]              

    > Output: list of 6mA modification from the combinations of motifs

    * [Included in modification_levels folder per barcode (wig format) for both plus and minus strand]


4) Data Analysis

    For the data analysis we use the outputs from modification_levels 
    Here in order to compare all possible motifs and barcodes combinations we have written some functions to compare them using the R scripting language. 

    For each of the comparison we cannot assume that they follow a normal distribution pattern so we
    decided to perform a non-parametric test with Wilcoxon also called Mann–Whitney U test (similar to a t-test where you compare mean distributions).
    We also do so for all the comparisons. We set a default parameter of two - sided where we do not specify which of the two distributions is greater than the other. 
    For example the idea is to see a significant difference with barcode05 with CYTANNNNNNVTGA against barcode03 with motif CYTANNNNNNVTGA. 
    Assuming that barcode05 did not show any 6mA for that motif. 

    * [***Important set working directory where the modification levels files are]
    Please open the script script.R in Rstudio in case you wanna re-run the test comparisons. 

    > Output: comparison file with all barcodes and all motifs

    * [Included in the main directory res.out]

