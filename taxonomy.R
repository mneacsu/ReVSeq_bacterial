pacman::p_load(rgbif)

genera <- read.csv('detected_genera.csv', header=T)

genera_unique <- unique(genera$genus)
  
genera_unique <- list_rbind(lapply(genera_unique, name_backbone))

genera_unique<-genera_unique[genera_unique$kingdom=='Bacteria','genus']

genera_unique<-unique(drop_na(genera_unique))

write.csv(genera_unique,'bacterial_genera.tsv', row.names = F, quote=F)
