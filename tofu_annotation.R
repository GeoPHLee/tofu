###annotation
library(biomaRt)
#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

annotation_dataframe=function(ase_set,parallel_core){
  parallel_core=as.numeric(parallel_core)
  
  if(parallel_core==1){
    gene_name_list=lapply(ase_set$event_name,gene_name_for_ase)
    gene_name_list=unlist(gene_name_list,recursive = FALSE)
    ase_set=cbind(ase_set,gene_name_list)
  }else{
    library(parallel)
    core_max=detectCores()-1
    if(core_max>=parallel_core){
      cl <- makeCluster(getOption("cl.cores", parallel_core))
      clusterExport(cl,"getBM")
      clusterExport(cl,"ensembl")
      clusterExport(cl,"position_praser_for_SE")
      clusterExport(cl,"position_praser_for_MXE")
      clusterExport(cl,"position_praser_for_RI")
      clusterExport(cl,"position_praser_for_A5SS")
      clusterExport(cl,"position_praser_for_A3SS")
      clusterExport(cl,"gene_name_for_ase")
      gene_name_list=parLapply(cl,ase_set$event_name,gene_name_for_ase)
      stopCluster(cl)
      gene_name_list=unlist(gene_name_list,recursive = FALSE)
      ase_set=cbind(ase_set,gene_name_list)
    }else{
      print("parallel_core required  is more than core_max,less!")
      return(NA)
    }
  }
  
  return(ase_set)
}
add_asetype_to_dataframe=function(ase_set){
  ase_type_list=lapply(ase_set$event_name,position_praser_for_all_type_return_ase_name)
  ase_type_list=unlist(ase_type_list,recursive = FALSE)
  ase_set=cbind(ase_set,ase_type_list)
  return(ase_set)
}

gene_name_for_ase=function(event_name_list){
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
    strand=-1
    #    # print(strand)
    #    #print("positive")
    #    strand=1
  }  
  if(strand=="+"){
    # print(strand)
    #print("positive")
    strand=1
  }
  
  chrosome_coordinate=list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  gene_name=getBM(attributes = c('hgnc_symbol'), filters = c('chromosome_name','start','end','strand'),values = chrosome_coordinate, mart = ensembl)
  #gene_name=getBM(attributes = c('mgi_symbol'), filters = c('chromosome_name','start','end','strand'),values = list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand), mart = ensembl)
  gene_name=unlist(gene_name)
  how_many_name=length(gene_name)
  if(how_many_name==1){
    return(gene_name)
  }else{
    if(how_many_name==0){
      gene_name=NA
      return(gene_name)
    }else{
      gene_name=Reduce(function(...) paste(...,sep=" "),gene_name)
      return(gene_name)
    }
  }
  
}
position_praser_for_A3SS=function(event_list){
  #A3SS_chr10:101458291:101458615:+@chr10:101460624|101460730:101460815:+
  all_position <- unlist(strsplit(unlist(as.vector(event_list)),split='@'))
  three_as_region=all_position[2]
  three_as_region=unlist(strsplit(unlist(as.vector(three_as_region)),split=':'))
  chr<- three_as_region[1]
  strand<- three_as_region[4]
  three_out=unlist(strsplit(three_as_region[2],split="\\|"))
  five_primer_exon_mRNA_start=three_out[1]
  three_primer_exon_end=three_out[2]
  if(strand =="+") mRNA_position<- list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  if(strand =="-") mRNA_position<- list(chr,three_primer_exon_end,five_primer_exon_mRNA_start,strand)
  return(mRNA_position)   
}
position_praser_for_A5SS=function(event_list){
  #A5SS_chr15:48623621:48623676|48623992:+@chr15:48624465:48624603:+
  all_position <- unlist(strsplit(unlist(as.vector(event_list)),split='@'))
  five_as_region=all_position[1]
  five_as_region=unlist(strsplit(unlist(as.vector(five_as_region)),split=':'))
  chr<- five_as_region[1]
  strand<- five_as_region[4]
  five_out=unlist(strsplit(five_as_region[3],split="\\|"))
  five_primer_exon_mRNA_start=five_out[1]
  three_primer_exon_end=five_out[2]
  if(strand =="+") mRNA_position<- list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  if(strand =="-") mRNA_position<- list(chr,three_primer_exon_end,five_primer_exon_mRNA_start,strand)
  return(mRNA_position)
  
}
position_praser_for_MXE=function(event_list){
  #MXE_chr11:43911275:43911378:+@chr11:43913591:43913679:+@chr11:43918745:43918904:+@chr11:43923066:43923275:+
  all_position <- unlist(strsplit(unlist(as.vector(event_list)),split='@'))
  region=all_position[2]
  region=unlist(strsplit(unlist(as.vector(region)),split=':'))
  chr<-region[1]
  strand<-region[4]
  five_primer_exon_mRNA_start<-region[2]
  three_primer_exon_end<-region[3]
  if(strand =="+") mRNA_position<- list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  if(strand =="-") mRNA_position<- list(chr,three_primer_exon_end,five_primer_exon_mRNA_start,strand)
  return(mRNA_position)
}
position_praser_for_RI=function(event_list){
  #RI_chr12:125609251-125609315:+@chr12:125609448-125609570:+

  all_position <- unlist(strsplit(unlist(as.vector(event_list)),split=':'))
  chr<-all_position[1]
  strand<-all_position[5]
  five_primer_exon_mRNA_start<-unlist(strsplit(unlist(as.vector(all_position[2])),split='-'))[2]
  three_primer_exon_end<-unlist(strsplit(unlist(as.vector(all_position[4])),split='-'))[1]
  if(strand =="+") mRNA_position<- list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  if(strand =="-") mRNA_position<- list(chr,three_primer_exon_end,five_primer_exon_mRNA_start,strand)
  return(mRNA_position)

  
}

position_praser_for_SE=function(event_list){
  #SE_chr4:106067842:106068136:+@chr4:106111517:106111643:+@chr4:106155054:106158508:+
  all_position <- unlist(strsplit(unlist(as.vector(event_list)),split='@'))

  SE_region=unlist(strsplit(unlist(as.vector(all_position[2])),split=':'))
  chr<-SE_region[1]
  strand<-SE_region[4]
  five_primer_exon_mRNA_start<-SE_region[2]
  three_primer_exon_end<-SE_region[3]
  if(strand =="+") mRNA_position<- list(chr,five_primer_exon_mRNA_start,three_primer_exon_end,strand)
  if(strand =="-") mRNA_position<- list(chr,three_primer_exon_end,five_primer_exon_mRNA_start,strand)
  return(mRNA_position)
  
  
  
}
position_praser_for_all_type_return_ase_name=function(event_name_list){
  ase <- unlist(strsplit(as.character(event_name_list),'\\_'))
  ase_type=ase[[1]]
  ase_event=ase[2]
  return(ase_type)
}
#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
cat("Example:\n\tase1_ase2=annotation_dataframe(ase1_ase2,cores)")





