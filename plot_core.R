ase_mitrix=select(ase_table,-1)

ase_mitrix2=mutate(ase_mitrix,sd=rowSds(as.matrix(ase_mitrix)))
ase_table2=ase_table[which(ase_mitrix2$sd>0.04),]
ase_40_m=ase_table2[-1]

library(gplots)
heatmap.2(as.matrix(ase_40_m),Rowv = FALSE,trace ="none" ,margins = c(10, 5),
offsetCol = 0,density.info = "none")

ggplot(melt(t(ase_40_m)), aes(x=Var2, y=Var1, fill=-log(value)))+geom_tile(color="white", size=0.01)+scale_fill_gradient(low='white', high='black')

library(ape)

ddata <- hclust(dist(t(as.matrix(ase_40_m))))
ddata <- hclust(dist(t(as.matrix(ase_40_m))),method="average")
ddata$labels=gsub("embryo","e",ddata$labels)
plot(as.phylo(ddata), type = "fan")


ase_table2=ase_table[which(ase_mitrix2$sd>0.02),]
ase_40_m=ase_table2[-1]
ddata <- hclust(dist(t(as.matrix(ase_40_m))),method="average")
gsub("embryo","e",ddata$labels)
ddata$labels=gsub("embryo","e",ddata$labels)
plot(as.phylo(ddata), type = "fan")

ggplot(melt(t(short_ase2[,2:37])), aes(Var2,Var1,fill=value,ylab="Early Embryo"))+geom_tile(color="white", size=0)+scale_fill_gradient(low='white', high='red',limits=seq(0,1),breaks=seq(0, 1,0.1))+scale_x_continuous("Event Index", breaks=seq(1, 8, 1))+theme(panel.grid =element_blank(),panel.border = element_blank())+theme_bw()+theme(axis.ticks = element_blank()) + theme(panel.grid =element_blank())+  theme(panel.border = element_blank())


short_ase=ase_table2_annotated[8:16,]

ase_table8=ase_table[which(ase_table$event_name%in%as.vector(short_ase$event_name)),]


