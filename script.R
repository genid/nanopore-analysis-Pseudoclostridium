
library(data.table)
setwd("/run/user/1005/gvfs/sftp:host=10.96.25.183/media/disk1/nanopore/Joyshree/project_2/guppy_basecaller/workspace/nanopore_pipeline/modification_levels/")
### Important set working directory where the modification levels files are
#setwd("/media/disk1/nanopore/Joyshree/project_2/guppy_basecaller/workspace/nanopore_pipeline/modification_levels/")
files = list.files(pattern = "*.out")

names = unlist(lapply(strsplit(files, split = "[.]"), function(x) paste(x[c(1, 3, 4)], collapse = "_")))

for(i in 1:length(files))
{
  assign(x = names[i], value = as.vector(read.table(files[i], header = F)$V1))
  print(i)
}

barcode_list <- lapply(1:length(names), function(x) get(names[x]))
names(barcode_list) <- names
head(barcode_list)

barcodes <- unlist(lapply(strsplit(names, split = "_"), function(x) x[1]))
motifs <- unlist(lapply(strsplit(names, split = "_"), function(x) x[4]))

make_comparison <- function(barcode1, barcode2, motif, barcodes, motifs, barcode_list)
{
  wilcox.test(x = unlist(barcode_list[(barcodes == barcode1) & (motifs == motif)]), 
              y = unlist(barcode_list[(barcodes == barcode2) & (motifs == motif)]), 
              alternative = "two.sided")$p.value
}

test_all_comparison <- function(barcodes, motifs, barcode_list)
{
  # Initialize list of matrices
  result_list <- list(length = length(unique(motifs)))
  for(i in 1:length(unique(motifs)))
  {
    result_list[[i]] <- matrix(nrow = length(unique(barcodes)), ncol = length(unique(barcodes)))
    rownames(result_list[[i]]) <- unique(barcodes)
    colnames(result_list[[i]]) <- unique(barcodes)
  }
  
  # Make comparisons
  for(i in 1:length(unique(motifs)))
  {
    for(j in 1:length(unique(barcodes)))
    {
      for(k in 1:length(unique(barcodes)))
      {
        tmp <- make_comparison(unique(barcodes)[j], unique(barcodes)[k], unique(motifs)[i], barcodes, motifs, barcode_list)
        result_list[[i]][j,k] <- tmp
      }
    }
    print(i)
  }
  names(result_list) <- unique(motifs)
  return(result_list)
}

res <- test_all_comparison(barcodes, motifs, barcode_list)
p.adjust(res$CACNNNNNNNTNGC, method = "bonferroni")
p.adjust(res$DGAGNNNNATC, method = "bonferroni")
res$GATNNNNCTC < 0.05
p.adjust(res$GATNNNNCTC, method = "bonferroni")
p.adjust(res$GATNNNNCTC, method = "BH")
p.adjust(res$GCNANNNNNNNGTG, method = "bonferroni")
p.adjust(res$TCABNNNNNNTARG, method = "bonferroni")

res$CACNNNNNNNTNGC <- p.adjust(res$CACNNNNNNNTNGC, method = "bonferroni")
res$DGAGNNNNATC <- p.adjust(res$DGAGNNNNATC, method = "bonferroni")
res$GATNNNNCTC <- p.adjust(res$GATNNNNCTC, method = "bonferroni")
res$GCNANNNNNNNGTG <- p.adjust(res$GCNANNNNNNNGTG, method = "bonferroni")
res$TCABNNNNNNTARG <- p.adjust(res$TCABNNNNNNTARG, method = "bonferroni")

p.vals.GATNNNNCTC <- c(1.100257e-21, 7.877415e-21,1.958053e-21,1.045101e-20,9.153713e-01,8.831269e-01,7.556155e-01,9.760171e-01,8.467131e-01,8.611344e-01)
p.vals.GATNNNNCTC + 1
p.adjust(p.vals.GATNNNNCTC, method = "bonferroni")
write.csv(file = "../GATNNNNCTC.corrected.txt",res$GATNNNNCTC, quote = F)
#write.csv(file = "res.out", x=res)

