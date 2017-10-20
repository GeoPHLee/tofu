inner_join2=function(x,y,by="event_name") inner_join(x,y, by="event_name")

one_sample_analysis=function(one_sample_name,inter_replicate_cut_off,replicate_times,name_of_sample){
  ##
  if(replicate_times==1){
    one_sample=one_replicate_all_ase(one_sample_name)
    return(one_sample)
  }else{
    one_replicate=map(one_sample_name,one_replicate_all_ase)
    one_sample  = reduce(one_replicate,inner_join2)
    one_sample_with_sd_filtered=filter_replicates(one_sample,cut_off,name_of_sample)
    return(one_sample_with_sd_filtered)
  }
  
  
}


one_replicate_all_ase=function(replicate_name){
  #EXAMPLE:all_ase_in_one_replicate=list("Late-blastocyst_embryo2_Cell1_SE_events.miso_summary","Late-blastocyst_embryo2_Cell1_RI_events.miso_summary","Late-blastocyst_embryo2_Cell1_A5SS_events.miso_summary","Late-blastocyst_embryo2_Cell1_A3SS_events.miso_summary","Late-blastocyst_embryo2_Cell1_MXE_events.miso_summary")
  ase_type<-c("_SE_events.miso_summary","_MXE_events.miso_summary","_A3SS_events.miso_summary","_A5SS_events.miso_summary","_RI_events.miso_summary")
  all_ase_in_one_replicate=paste(replicate_name,ase_type,sep = "")
  #print(all_ase_in_one_replicate)
  #print("OK")
  ase_sets=map_dfr(all_ase_in_one_replicate,summary_file_process_for_one_ase)
  ##one_rep_all_ase <- bind_rows(ase_sets)
  return(ase_sets)
}

summary_file_process_for_one_ase=function(sample_replicate){
  
  #EXAMPLE:all_replicates="Late-blastocyst_embryo2_Cell1_SE_events.miso_summary"
  ##get sample information
      #file_name  = as.character(sample_replicate)
  file_name=sample_replicate
      #which sample i am working onsy
      #print(file_name)
      #print("ok1")
  ##read table
  psi_table = read_tsv(file_name,col_names = TRUE,progress = FALSE)
  ##select all psi or according to a given cut_off <= 1.0
  psi_table = psi_table[which(abs(psi_table$ci_high-psi_table$ci_low)<=delta_cutoff),]

      #print("read ok")
  ##select all psi or according to a given cut_off <= 1.0
      #if (delta_cutoff == 1.0){
      #  tem_result = matrix(tem_result,nc=3,byrow = T)
      #}else{
      #  tem_result= t(as.data.frame(lapply(as.matrix(tem_result),unlist)))
      # }
  psi_table = select(psi_table, event_name,miso_posterior_mean)
  ##get metadata
  orig_table_name  =  substr(file_name,1,stop = as.numeric(gregexpr("\\_events\\.miso_summary",file_name))-1)
  ase_event_type=substr(file_name,regexpr("[A-Z|a-z|0-9]+\\_events",file_name),regexpr("events",file_name)-2)
  ##mark ase type
  psi_table=mutate(psi_table,event_name=paste(ase_event_type,event_name,sep = "_"))
  ##mark sample name
  sample_name=substr(orig_table_name,1,stop=as.numeric(gregexpr(ase_event_type,orig_table_name))-2)
  names(psi_table)[2]=sample_name
      # psi_list = unlist(tem_result[,2])
      #  event_name_list = unlist(tem_result[,1])
      #    delta_ci_list =  as.numeric(unlist(tem_result[,3]))
      #psi_flag ="psi_"
      #  delta_flag = "delta_ci_"
      #  psi_table = table_process(event_name_list,psi_list,orig_table_name,psi_flag,ase_event_type)
      #    delta_psi_table =table_process(event_name_list,delta_ci_list,orig_table_name,delta_flag)
      #    result = list(psi=psi_table,delta_psi=delta_psi_table)
 return (psi_table)
}
