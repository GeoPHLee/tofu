sample_list<-"New_Sample_List.csv"
sample_list_file <- read_csv(sample_list,col_names = TRUE)

group_list=unique(sample_list_file$Group)


all_sample_ase_list<-parApply(cl,sample_list_file,1,one_sample_analysis)
ase_table=purrr::reduce(all_sample_ase_list,inner_join2)
ase_matrix=ase_table[,-1]
ase_dist=dist(t(as.matrix(ase_matrix)))
tag_genator=function(sample_info){
  tag=rep(paste(sample_info[[1]],sample_info[[2]],sep="_"),sample_info[[3]])
  return(tag)
}
group_level_rep=apply(sample_list_file,1,tag_genator)

ge_result=pca_for_dataframe(ase_table2,2,unlist(group_level_rep))
#group_list=c("Oocyte" ,"Zygote","2-cell","4-cell","8-cell","Morulae","Late-blastocyst")
#group_list=c("oocyte" ,"zygote","pronuclei","2-cell")
group_list=c("Oocyte" ,"Zygote","2-cell","4-cell","8-cell","Morulae","Late-blastocyst")
sample_list<-"sample_list.csv"
sample_list_file <- read.csv(sample_list,header = TRUE,stringsAsFactors = FALSE)
each_sample_replicate<-lapply(apply(sample_list_file,1,as.list),unlist)
source("tofu_misoProcessing.R")

all_sample_ase_list<-analysis_for_all_ase(each_sample_replicate,group_list)
ase_table=Reduce(function(...) merge(..., by="event_name_list"),unlist(all_sample_ase_list,recursive = FALSE))

source("tofu_rsemProcessing.R")
ge_summary_exits=FASLE
ge_summary_file=""
if(ge_summary){
	gene_expression_table=read.csv(ge_summary_file)
}else{
  all_tj_ge_list<-parApply(cl,sample_list_file,1,one_sample_analysis)
	ge_table=purrr::reduce(all_tj_ge_list,inner_join2)
}
#Define group_level
##group_level=c("oocyte","oocyte","oocyte","zygote","zygote","zygote","2-cell-cell1","2-cell-cell1","2-cell-cell1","2-cell-cell2","2-cell-cell2","2-cell-cell2")
##group_level=factor(group_level,ordered = TRUE,levels = c("oocyte","zygote","2-cell-cell1","2-cell-cell2"))
group_level=sample_list_file[1,]
group_level=factor(group_level,ordered = TRUE)
###gobal level:PCA 
source("tofu_PCA.R")

###local level
source("tofu_groupCompare.R")
oocyte=oztf_all_list[[1]]
zygote=oztf_all_list[[2]]
two_cell=oztf_all_list[[3]]
four_cell=oztf_all_list[[4]]
eight_cell=oztf_all_list[[5]]
morulae=oztf_all_list[[6]]
lb=oztf_all_list[[7]]

oocyte1_zygote=allg1_allg2(oocyte1,zygote,0.2)
zygote_two_cell=allg1_allg2(zygote,two_cell,0.2)
two_cell_four_cell=allg1_allg2(two_cell,four_cell,0.2)
four_cell_eight_cell=allg1_allg2(four_cell,eight_cell,0.2)
eight_cell_morulae=allg1_allg2(eight_cell,morulae,0.2)
morulae_lb=allg1_allg2(morulae,lb,0.2)

###annotation and detail analysis
source("tofu_annotation_final.R")
source("tofu_GEvsASE.R")
