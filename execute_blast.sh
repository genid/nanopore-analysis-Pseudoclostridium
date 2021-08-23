blastn -query motifs/CACNNNNNNNTNGC.txt -db ecoli/HST04 -out motifs/CACNNNNNNNTNGC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 20;
blastn -query motifs/DGAGNNNNATC.txt -db ecoli/HST04 -out motifs/DGAGNNNNATC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 20;
blastn -query motifs/GCNANNNNNNNGTG.txt -db ecoli/HST04 -out motifs/GCNANNNNNNNGTG.out -task blastn-short -outfmt 7 -perc_identity 100 -num_threads 20;
blastn -query motifs/TCABNNNNNNTARG.txt -db ecoli/HST04 -out motifs/TCABNNNNNNTARG.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 20;
blastn -query motifs/CYTANNNNNNVTGA.txt -db ecoli/HST04 -out motifs/CYTANNNNNNVTGA.out -task blastn-short -outfmt 7 -perc_identity 100 -num_threads 20;
blastn -query motifs/GATNNNNCTC.txt -db ecoli/HST04 -out motifs/GATNNNNCTC.out -task blastn-short -outfmt 7 -perc_identity 100 -evalue 7 -num_threads 20;
