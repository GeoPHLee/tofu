#environment variable 
##need_replicate(TRUE/FALSE)
inner_join2=function(x,y,by="gene_id") inner_join(x,y, by="gene_id")
full_join2=function(x,y,by="gene_id") full_join(x,y, by="gene_id")
library(parallel)
library(tidyverse)
library(matrixStats)
inter_replicate_cut_off=0.05
need_replicate=TRUE
delta_cutoff=0.3
lowest_TPM=5
sd_cutoff=1
need_filter=FALSE
repli_name="_embryo"
cl <- makeCluster(getOption("cl.cores", parallel_core))
#clusterExport(cl,"tidyverse")
clusterExport(cl,"inner_join2")
clusterExport(cl,"full_join2")

#clusterExport(cl,"matrixStats")
clusterExport(cl,"repli_name")
clusterExport(cl,"need_replicate")
clusterExport(cl,"delta_cutoff")
clusterExport(cl,"lowest_TPM")
clusterExport(cl,"sd_cutoff")
clusterExport(cl,"need_filter")
#delta_cutoff=1
clusterExport(cl,"inter_replicate_cut_off")
#inter_replicate_cut_off=0.05

clusterExport(cl,"one_sample_analysis")
clusterExport(cl,"filter_replicates")
clusterExport(cl,"resm_table_file_process_for_one_replicate")
clusterExport(cl,"one_replicate_ge")
#clusterExport(cl,"one_replicate_ge")

one_sample_analysis=function(group_condition_replicate_times){
  ##one_sample_name=c("Late-blastocyst_embryo1_Cell1","Late-blastocyst_embryo2_Cell1","Late-blastocyst_embryo3_Cell1")
  library(tidyverse)
  group=group_condition_replicate_times[[1]]
  condition=group_condition_replicate_times[[2]]
  replicate_times=group_condition_replicate_times[[3]]
  
  
  
  replicate_id<-seq(1:replicate_times)
  all_replicate_name<-paste(repli_name,replicate_id,"_",sep="")
  one_sample_name=paste(group,all_replicate_name,condition,sep="")
  if(replicate_times==1){
    one_sample=one_replicate_ge(one_sample_name)
    return(one_sample)
  }else{
    all_replicates=purrr::map(one_sample_name,one_replicate_ge)
    one_sample=purrr::reduce(all_replicates,full_join2)
  #  one_sample=purrr::reduce(all_replicates,inner_join2)
#one_sample=do.call(inner_join2,all_replicates)
    if(need_filter==TRUE){
        one_sample_with_sd_filtered=filter_replicates(one_sample,inter_replicate_cut_off,one_sample_name)
        return(one_sample_with_sd_filtered)
        }else{
        return(one_sample)
        }
  }
}

filter_replicates=function(sample_results,sd_cutoff,name_of_sample){
  ##assign environment varaible  need_replicate
  library(matrixStats)
  item_list<-dplyr::select(sample_results,-1)
  one_sample=mutate(sample_results,psi_mean=round(rowMeans(item_list),digits = 2),intersample_sd=rowSds(as.matrix(item_list)/round(rowMeans(item_list),digits = 2), na.rm=TRUE))
  sample_result_mean_sd_filtered<-one_sample[which(one_sample$intersample_sd<sd_cutoff),]
  
  #  print(name_of_sample)
  #  print(names(sample_result_mean)[2])
  data_frame_length=length(sample_result_mean_sd_filtered)
  if(need_replicate==FALSE){
    names(sample_result_mean_sd_filtered)[data_frame_length-1]<-name_of_sample
    sample_result_mean_sd_filtered=sample_result_mean_sd_filtered[,c(1,data_frame_length-1)]
  }else{
    sample_result_mean_sd_filtered=sample_result_mean_sd_filtered[1:(data_frame_length-2)]
  }
  return(sample_result_mean_sd_filtered)
}
one_replicate_ge=function(replicate_name){
  ##replicate_name="Late-blastocyst_embryo2_Cell1"
  #EXAMPLE:all_ase_in_one_replicate=list("Late-blastocyst_embryo2_Cell1_SE_events.miso_summary","Late-blastocyst_embryo2_Cell1_RI_events.miso_summary","Late-blastocyst_embryo2_Cell1_A5SS_events.miso_summary","Late-blastocyst_embryo2_Cell1_A3SS_events.miso_summary","Late-blastocyst_embryo2_Cell1_MXE_events.miso_summary")
  ge_name<-"_genecode.genes.results"
  ge_file_name=paste(replicate_name,ge_name,sep = "")
  #print(all_ase_in_one_replicate)
  #print("OK")
  ge_sets=resm_table_file_process_for_one_replicate(ge_file_name)
  ##one_rep_all_ase <- bind_rows(ase_sets)
  return(ge_sets)
#    ase_type<-c("#_SE_events.miso_summary","#_MXE_events.miso_summary","#_A3SS_events.miso_summary",#"_A5SS_events.miso_summary"#,"_RI_events.miso_summary")
#  all_ase_in_one_replicate=paste(#replicate_name,ase_type,sep =# "")
#  #print(all_ase_in_one_replicate)
#  #print("OK")
#  ase_sets=map_dfr(#all_ase_in_one_replicate,summ#ary_file_process_for_one_ase)
#  ##one_rep_all_ase <- bind_rows(#ase_sets)
#  return(ase_sets)
}


resm_table_file_process_for_one_replicate=function(sample_replicate){

    #get sample information
    file_name  = sample_replicate
    orig_table_name  =  substr(file_name,1,stop = as.numeric(gregexpr("genecode\\.genes\\.results",file_name))-1)
    #which sample i am working on
    print(file_name)
    #print("ok1")
    #read table
    resm_table = read_tsv(file_name,col_names = TRUE,progress = FALSE)
    resm_table = resm_table[resm_table$TPM>lowest_TPM,]
    resm_table = dplyr::select(resm_table, gene_id,TPM)
    print("read ok")
    flag ="TPM_"
    list_newname = paste(flag,orig_table_name,sep="")
    names(resm_table)[2]=list_newname
    return (resm_table)
}
