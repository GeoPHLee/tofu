library(pca3d)
data(metabo)
pca <- prcomp(metabo[,-1], scale.=TRUE)
gr <- factor(metabo[,1])
summary(gr)

four_cell
two_cell1_ase_name=all_sample_ase_list[[33]][,1]
two_cell1_ase_name_list=as.list(unlist(two_cell1_ase_name))
two_cell2_ase_name=all_sample_ase_list[[34]][,1]
two_cell2_ase_name_list=as.list(unlist(two_cell2_ase_name))

all_2_cell_ase=purrr::reduce(list(all_sample_ase_list[[33]],all_sample_ase_list[[34]]),full_join)
all_2_cell_ase=do.call(full_join,list(all_sample_ase_list[[33]],all_sample_ase_list[[34]]))

twocell_1ase_unique=all_two_cell_ase[which(is.na(all_two_cell_ase$`2-cell_embryo1_Cell2`)==TRUE),]
twocell_2ase_unique=all_two_cell_ase[which(is.na(all_two_cell_ase$`2-cell_embryo1_Cell1`)==TRUE),]

gene2cell1 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(twocell_1ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db", IDs2Add = "symbol")
gene2cell2 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(twocell_2ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db", IDs2Add = "symbol")
venn.diagram(x=list(`#1-2_cell`=na.omit(gene2cell1$symbol),`#2-2_cell`=na.omit(gene2cell2$symbol)),filename="ase_diff_gene2cell.png",fill=c(colors()[148], colors()[589]),imagetype="png")

two1gu=base::setdiff(na.omit(gene2cell1$symbol),na.omit(gene2cell2$symbol))
two2gu=base::setdiff(na.omit(gene2cell2$symbol),na.omit(gene2cell1$symbol))
write_lines(na.omit(gene2cell1$symbol),"2cell_1ase_unique_gene_name")
write_lines(na.omit(gene2cell2$symbol),"2cell_2ase_unique_gene_name")


venn.diagram(x=list(`#1-4_cell`=na.omit(gene4cell1$symbol),`#2-4_cell`=na.omit(gene4cell2$symbol),`#3-4_cell`=na.omit(gene4cell3$symbol),`#4-4_cell`=na.omit(gene4cell4$symbol)),filename="ase_in_gene_diff4cell.png",fill=c("cornflowerblue","green","yellow","darkorchid1"),imagetype="png")


four_cell1_ase_name=all_four_list[[1]][,1]
four_cell1_ase_name_list=as.list(unlist(four_cell1_ase_name))

four_cell2_ase_name=all_four_list[[2]][,1]
four_cell2_ase_name_list=as.list(unlist(four_cell2_ase_name))

four_cell3_ase_name=all_four_list[[3]][,1]
four_cell3_ase_name_list=as.list(unlist(four_cell3_ase_name))
four_cell4_ase_name=all_four_list[[4]][,1]
four_cell4_ase_name_list=as.list(unlist(four_cell4_ase_name))

venn.diagram(x=list(`#1-4_cell`=four_cell1_ase_name_list,`#2-4_cell`=four_cell2_ase_name_list,`#3-4_cell`=four_cell3_ase_name_list,`#4-4_cell`=four_cell4_ase_name_list),filename="ase_diff4cell.png",fill=c("cornflowerblue","green","yellow","darkorchid1"),imagetype="png")

venn.diagram(x=list(`#1-2_cell`=two_cell1_ase_name_list,`#2-2_cell`=two_cell2_ase_name_list),filename="ase_diff2cell.png",fill=c(colors()[148], colors()[589]),imagetype="png")



all_four_cell_ase=purrr::reduce(all_four_list,full_join)
fourcell_1ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)==TRUE),]










fourcell_1ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)==TRUE),]
View(fourcell_1ase_unique)
fourcell_2ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell1`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)==TRUE),]
fourcell_3ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell1`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)==TRUE),]
fourcell_4ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell1`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)==TRUE),]
View(fourcell_1ase_unique)
gene4cell1 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_1ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_1ase_unique_gene_name")
gene4cell2 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_2ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_2ase_unique_gene_name")
gene4cell3 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_3ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_3ase_unique_gene_name")
gene4cell4 <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_4ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_4ase_unique_gene_name")
View(all_four_cell_ase)


venn.diagram(x=list(`#1-4_cell`=na.omit(gene4cell1$symbol),`#2-4_cell`=na.omit(gene4cell2$symbol),`#3-4_cell`=na.omit(gene4cell3$symbol),`#4-4_cell`=na.omit(gene4cell4$symbol)),filename="ase_in_gene_diff4cell.png",fill=c("cornflowerblue","green","yellow","darkorchid1"),imagetype="png")


fourcell_ase_inner=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell1`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)==FALSE),]
overlaps.anno <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_ase_inner$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_ase_inner")
overlaps.anno <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_ase_inner$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(unique(na.omit(overlaps.anno$symbol)),"fourcell_ase_inner")



four1g=as.list(unique(na.omit(gene4cell1$symbol)))
four2g=as.list(unique(na.omit(gene4cell2$symbol)))
four3g=as.list(unique(na.omit(gene4cell3$symbol)))
four4g=as.list(unique(na.omit(gene4cell4$symbol)))


four1gu=base::setdiff(four1g,purrr::reduce(list(four2g,four3g,four4g),union))
four2gu=base::setdiff(four2g,purrr::reduce(list(four1g,four3g,four4g),union))
four3gu=base::setdiff(four3g,purrr::reduce(list(four2g,four1g,four4g),union))
four4gu=base::setdiff(four4g,purrr::reduce(list(four2g,four3g,four1g),union))



cor_to_range=function(event_name_list){
library(GenomicRanges)
ase <- unlist(strsplit(as.character(event_name_list),'\\$'))
ase_type=ase[[1]]
ase_event=ase[2]
mRNA_position=switch(ase_type,
SE =position_praser_for_SE(ase_event),
A5SS=position_praser_for_A5SS(ase_event),
A3SS=position_praser_for_A3SS(ase_event),
MXE=position_praser_for_MXE(ase_event),
RI=position_praser_for_RI(ase_event)
)
print(mRNA_position)
chr <-  unlist(mRNA_position[1])
chr<-substr(chr,4,nchar(chr)+1)
if(chr=="Y"|chr=="X"){
chr=chr
}else{
chr<-as.numeric(chr)
}
five_primer_exon_mRNA_start <- as.numeric(as.character(mRNA_position[2]))
three_primer_exon_end <- as.numeric(as.character(mRNA_position[3]))
strand <- as.character(mRNA_position[4])
#  if(strand=="+"){
#    # print(strand)
#    #print("positive")
#    strand=1
#  }else{
#    #  print("no")
#  }
#  if(strand=="-"){
#    #print(strand)
#    # print("minus")
#    strand=-1
#  }
# chrosome_coordinate=list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
#  gene_name=getBM(attributes = c('hgnc_symbol'), filters = c('chromosome_name','start','end','strand'),values = chrosome_coordinate, mart = ensembl)
#gene_name=getBM(attributes = c('mgi_symbol'), filters = c('chromosome_name','start','end','strand'),values = list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand), mart = ensembl)
query_range <- GRanges(chr,IRanges(five_primer_exon_mRNA_start,three_primer_exon_end),strand)
return(query_range)
}