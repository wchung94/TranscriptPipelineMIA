#""" Create non-relative and relative sized venndiagrams
#WY Chung
#"""

if(length(.libPaths()) == 1){
    # We're in Rscript.exe
    possible_lib_paths <- file.path(Sys.getenv(c('USERPROFILE','R_USER')),
                                    "R","~/R/x86_64-pc-linux-gnu-library/3.3",
                                    paste(R.version$major,
                                             substr(R.version$minor,1,1),
                                             sep='.'))
    indx <- which(file.exists(possible_lib_paths))
    if(length(indx)){
       .libPaths(possible_lib_paths[indx[1]])
    }
    # CLEAN UP
    rm(indx,possible_lib_paths)
}


#install.packages('VennDiagram')
library(VennDiagram)
#install.packages('eulerr')
library(eulerr)

mia_file = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast/ia_lines_exp.txt'
mia_lines_exp = read.table(mia_file, header=FALSE)


num_genes = nrow(mia_lines_exp)
CK_MJ = sum(mia_lines_exp$SigCK_MJ == 'no' )
CK_ET = sum(mia_lines_exp$SigCK_ET == 'no')
MJ_ET = sum(mia_lines_exp$SigMJ_ET == 'no' )
CK_MJ_ET = sum(mia_lines_exp$SigCK_ET == 'no' & mia_lines_exp$SigCK_MJ == 'no' & mia_lines_exp$SigMJ_ET == 'no')


# Plot venn diagram with numbers but without relative size
grid.newpage()
draw.triple.venn(area1 = num_genes, area2 = num_genes, area3 = CK_MJ + CK_ET, n12 = MJ_ET, n23 = CK_ET, n13 = CK_MJ, 
                 n123 = CK_MJ_ET, category = c("MJ", "ET", "CK"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))


# Plot venn diagram without numbers but with relative size
grid.newpage()
VennDiag_MJ <- euler(c("MJ" = num_genes-CK_MJ,  "CK" = num_genes-CK_MJ,   "MJ&CK" = CK_MJ ))
plot(VennDiag_MJ, counts = TRUE, print.mode=("raw"), 
     main = 'MJ vs CK', font=2, cex=1, alpha=0.5, lty = 'blank', category = c("MJ", "CK"),
     fill=c("skyblue", "mediumorchid"))

VennDiag_ET <- euler(c("ET" = num_genes-CK_ET,  "CK" = num_genes-CK_ET,   "ET&CK" = CK_ET ))
plot(VennDiag_ET, counts = TRUE, print.mode=("raw"), 
     main = 'ET vs CK', font=2, cex=1, alpha=0.5, lty = 'blank', category = c("ET", "CK"),
     fill=c("pink1", "mediumorchid"))

