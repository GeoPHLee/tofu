library(pca3d)
data(metabo)
pca <- prcomp(metabo[,-1], scale.=TRUE)
gr <- factor(metabo[,1])
summary(gr)

four_cell
two_cell1_ase_name=all_two_list[[2]][,1]
two_cell1_ase_name_list=as.list(unlist(two_cell1_ase_name))
two_cell2_ase_name=all_two_list[[2]][,1]
two_cell2_ase_name_list=as.list(unlist(two_cell2_ase_name))

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


twocell_1ase_unique=all_two_cell_ase[which(is.na(all_two_cell_ase$`2-cell_embryo1_Cell2`)==TRUE),]
twocell_2ase_unique=all_two_cell_ase[which(is.na(all_two_cell_ase$`2-cell_embryo1_Cell1`)==TRUE),]

all_four_cell_ase=purrr::reduce(all_four_list,full_join)
fourcell_1ase_unique=all_four_cell_ase[which(is.na(all_four_cell_ase$`4-cell_embryo1_Cell2`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell3`)&is.na(all_four_cell_ase$`4-cell_embryo1_Cell4`)==TRUE),]

twocell_2ase_unique=all_two_cell_ase[which(is.na(all_two_cell_ase$`2-cell_embryo1_Cell1`)==TRUE),]