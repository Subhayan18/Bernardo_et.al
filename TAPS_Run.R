# Collect the data set
lapply(c('tidyr','dplyr','data.table'),require, character.only=T)
file.name<-list.files()[grep(".table",list.files())]
df<-fread(file.name)

df$Chr<-paste0('chr',df$Chr)

probes <- df %>% select(Chromosome = Chr,
                     Start = Position,
                     End = Position,
                     Value = colnames(df)[grep('LRR',colnames(df))]) %>% 
  mutate(End = End + 50)
write.table(probes, "probes.txt", row.names=FALSE, quote=FALSE, sep='\t')

snps <- df %>% select(Chromosome = Chr,
              Start = Position,
              End = Position,
              Value = colnames(df)[grep('BAF',colnames(df))])
write.table(snps, "snps.txt", row.names=FALSE, quote=FALSE, sep='\t')

options(bitmapType='cairo')
TAPS::TAPS_plot()
