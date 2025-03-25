lapply(c('tidyr','dplyr'),require, character.only=T)

load("D:/BMC/Box/Helping with analysis/Carina/Data/PDX_CNV.RData")

all.names<-c(colnames(cnv_b1),colnames(cnv_b2),colnames(cnv_b3))
common.Name<-intersect(intersect(cnv_b1$Name,cnv_b2$Name),cnv_b3$Name)

cnv_b1_filtered <- cnv_b1 %>% 
  filter(Name %in% common.Name == TRUE) %>%  
  filter(Chr != 0) %>% 
  filter(Position != 0) %>% 
  mutate(idn = paste0(Chr,":",Position)) %>% 
  filter(duplicated(idn) == FALSE)

cnv_b2_filtered <- cnv_b2 %>% 
  filter(Name %in% common.Name == TRUE) %>%  
  filter(Chr != 0) %>% 
  filter(Position != 0) %>% 
  mutate(idn = paste0(Chr,":",Position)) %>% 
  filter(duplicated(idn) == FALSE)

cnv_b3_filtered <- cnv_b3 %>% 
  filter(Name %in% common.Name == TRUE) %>%  
  filter(Chr != 0) %>% 
  filter(Position != 0) %>% 
  mutate(idn = paste0(Chr,":",Position)) %>% 
  filter(duplicated(idn) == FALSE)

rm(list=c('cnv_b1','cnv_b2','cnv_b3'))
gc()

common.idn<-intersect(intersect(cnv_b1_filtered$idn,cnv_b2_filtered$idn),cnv_b3_filtered$idn)

cnv_b1_filtered <- cnv_b1_filtered %>%
  filter(idn %in% common.idn) %>% 
  arrange(Chr,Position)
cnv_b2_filtered <- cnv_b2_filtered %>%
  filter(idn %in% common.idn)%>% 
  arrange(Chr,Position)
cnv_b3_filtered <- cnv_b3_filtered %>%
  filter(idn %in% common.idn)%>% 
  arrange(Chr,Position)

#Now all variants are aligned, sorted across three data sets ((~700k variants dropeed out altogether)

cnv_filtered <-cbind (cnv_b1_filtered[,-c(1,ncol(cnv_b1_filtered))],
                    cnv_b2_filtered[,-c(1:4,ncol(cnv_b2_filtered))],
                    cnv_b3_filtered[,-c(1:4,ncol(cnv_b3_filtered))])

model.1 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"1")]))
model.2 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"2")]))
model.3 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"3")]))
model.4 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"4")]))
model.5 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"5")]))
model.6 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"6")]))
model.7 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"7")]))
model.8 <- cnv_filtered %>% select(c('Name','Chr','Position',all.names[startsWith(all.names,"8")]))

#Writing the model files in the server directories
write.table(model.1, "~/Carina_mice/model1/model1.table", row.names=FALSE, quote=FALSE)
write.table(model.2, "~/Carina_mice/model2/model2.table", row.names=FALSE, quote=FALSE)
write.table(model.3, "~/Carina_mice/model3/model3.table", row.names=FALSE, quote=FALSE)
write.table(model.4, "~/Carina_mice/model4/model4.table", row.names=FALSE, quote=FALSE)
write.table(model.5, "~/Carina_mice/model5/model5.table", row.names=FALSE, quote=FALSE)
write.table(model.6, "~/Carina_mice/model6/model6.table", row.names=FALSE, quote=FALSE)
write.table(model.7, "~/Carina_mice/model7/model7.table", row.names=FALSE, quote=FALSE)
write.table(model.8, "~/Carina_mice/model8/model8.table", row.names=FALSE, quote=FALSE)

