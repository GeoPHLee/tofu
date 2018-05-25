library(EnsDb.Hsapiens.v75) ##(hg19)
library(ChIPpeakAnno)
library(GenomicRanges)
cl <- makeCluster(getOption("cl.cores", parallel_core))
clusterExport(cl,"position_praser_for_SE")
clusterExport(cl,"position_praser_for_MXE")
clusterExport(cl,"position_praser_for_RI")
clusterExport(cl,"position_praser_for_A5SS")
clusterExport(cl,"position_praser_for_A3SS")
cor_to_range<-function(event_name_list){
  library(GenomicRanges)
  ase <- unlist(strsplit(as.character(event_name_list),'\\$'))
  ase_type=ase[[1]]
  print(ase)
  ase_event=ase[2]
  mRNA_position=switch(ase_type,
                       SE =position_praser_for_SE(ase_event),
                       A5SS=position_praser_for_A5SS(ase_event),
                       A3SS=position_praser_for_A3SS(ase_event),
                       MXE=position_praser_for_MXE(ase_event),
                       RI=position_praser_for_RI(ase_event)
  )
  print(mRNA_position)
  # chr <-  unlist(mRNA_position[1])
  chr<-substr(mRNA_position[1],4,nchar(mRNA_position[1])+1)
  if(grepl("^[[:digit:]]+$",chr)){
    chr<-as.numeric(chr)
  }else{
    chr<-chr
  }
  five_primer_exon_mRNA_start <- as.numeric(as.character(mRNA_position[2]))
  three_primer_exon_end <- as.numeric(as.character(mRNA_position[3]))
  strand <- as.character(mRNA_position[4])
  if(strand=="-"){
    temp=five_primer_exon_mRNA_start
    five_primer_exon_mRNA_start=three_primer_exon_end
    three_primer_exon_end=temp
    #    # print(strand)
    #    #print("positive")
    #    strand=1
  }
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

overlaps.anno <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,ase_tablpe$event_name,1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db", IDs2Add = "symbol")

write_lines(na.omit(overlaps.anno$symbol),"fourcell_1ase_unique_gene_name")



