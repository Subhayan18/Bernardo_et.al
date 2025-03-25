install.packages("TAPS", repos="http://R-Forge.R-project.org")
BiocManager::install("DNAcopy")
BiocManager::install("affxparser")
install.packages("fields")
install.packages("foreach")
install.packages("jpeg")
install.packages("xlsx")
require(dplyr)
require(readr)
require(TAPS)
model4_cnv <- read_csv("Data/model4_cnv.csv", 
                       col_types = cols(...1 = col_skip(), Index = col_integer()))
View(model4_cnv)


model4_cnv$Chr<-paste0("chr",model4_cnv$Chr)

S_pr =  model4_cnv %>% select('Chromosome' = Chr, 
                        'Start' = Position, 
                        'End' = Position, 
                        'Value' = `4P2C_1E.LRR`)
S_pr$End = S_pr$End+10

S_snp =  model4_cnv %>% select('Chromosome' = Chr, 
                                 'Start' = Position, 
                                 'End' = Position, 
                                 'Value' = `4P2C_1E.BAF`)
S_snp$End = S_snp$End+10

setwd("D:/BMC/Box/Helping with analysis/Carina/Model4D_example/4P2C1E/")
write.table(S_pr, "probes.txt",sep='\t',quote = F, row.names = F)
write.table(S_snp, "snps.txt",sep='\t',quote = F, row.names = F)
TAPS_plot()
