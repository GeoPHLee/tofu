overlaps.anno <- addGeneIDs(annotatePeakInBatch(fourlist,AnnotationData=annoData, output="overlapping") ,
                            "org.Hs.eg.db",
                            IDs2Add = "symbol")



cl <- makeCluster(getOption("cl.cores", parallel_core))
clusterExport(cl,"GenomicRanges")
clusterExport(cl,"position_praser_for_SE")
clusterExport(cl,"position_praser_for_MXE")
clusterExport(cl,"position_praser_for_RI")
clusterExport(cl,"position_praser_for_A5SS")
clusterExport(cl,"position_praser_for_A3SS")
test_four_ase=parApply(cl,four_cell1_ase_name,1, cor_to_range)


all_four_ase2=do.call(c,parApply(cl,four_cell1_ase_name,1, cor_to_range))
overlaps.anno <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,four_cell1_ase_name,1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db", IDs2Add = "symbol")

#all_four_ase2=Reduce(c,all_four_ase)
#all_four_ase2=purrr::reduce(all_four_ase,c)


all_four_ase2=do.call(c,all_four_ase)



overlaps.anno <- addGeneIDs(annotatePeakInBatch(do.call(c,parApply(cl,as.data.frame(fourcell_1ase_unique$event_name),1, cor_to_range)),AnnotationData=annoData, output="overlapping") ,"org.Hs.eg.db",                            IDs2Add = "symbol")
write_lines(na.omit(overlaps.anno$symbol),"fourcell_1ase_unique_gene_name")