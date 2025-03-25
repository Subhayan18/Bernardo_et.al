lapply(c('tidyr','dplyr','data.table'),require, character.only=T)
file.name<-list.files()[grep("model",list.files())]
df<-fread(file.name)
samples <- colnames(df)[grep("LRR",colnames(df))]
samples <- gsub("\\..*","",samples)

for (i in 1 : length(samples)) {
  assign(samples[i], df %>% select(c('Name','Chr','Position',colnames(df)[grep(samples[i],colnames(df))][1:3])))
  dir.create(samples[i])
  write.table(get(samples[i]), paste0(getwd(),"/",samples[i],"/",samples[i],".table"), row.names=FALSE, quote=FALSE)
}