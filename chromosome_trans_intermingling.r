## Find genome loci ( here beads in genome structural models) that are in close proximity to each other  
## Parameters: (1) PATH_nuc is the folder path containing files (hdf5 format) of calculated genome structures
## (2) proximity threshold can normally be set between 1.5 and 3.0; a threshold of 2 is mainly used in this study, 
##     but threshold as low as 1.2 or as high as 3.5 all gave similar results in terms of intermingling score analysis through different time points and the comparison between A/B compartment
##     It has the unit of bead radii in genome structure polymer models. 
Dis3D_proximal=function(PATH_nuc,PATH_output,proximity_threshold,model_num_total){
  library('rhdf5')
  q<-dir(PATH_nuc)
  k=1
  fraction_of_bead_pairs_consistent_within_distance_thresh=c()
  
  #####
  distance_function=function(vector1,vector2){
    sqrt((vector1[1]-vector2[1])^2+(vector1[2]-vector2[2])^2+(vector1[3]-vector2[3])^2)
  } 
  
  distdex<-function(i,j,n) #given row, column, and n, return index
    n*(i-1) - i*(i-1)/2 + j-i
  
  rowcol<-function(ix,n) { #given index, return row and column
    nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
    nc=n-(2*n-nr+1)*nr/2+ix+nr
    cbind(nr,nc)
  }
  #####
  
  for (file in q){
    
    chr_bead_num_each_file=c()
    
    a=h5read(paste(PATH_nuc,file,sep="/"), "/")
    for(model_num in 1:model_num_total){
      posi_matrix_each_model=c()
      for (chr_num in c(1:20)){
        if(chr_num %in% 1:19 ==TRUE ){ chr_num_old_label=as.character(chr_num)}
        if(chr_num ==20 ){ chr_num_old_label="X"}
        coor_chr_each_model=a[["structures"]][[1]][["coords"]][[chr_num_old_label]][,,model_num]
        if(model_num==1){chr_bead_num_each_file=c(chr_bead_num_each_file,dim(coor_chr_each_model)[2])}
        posi_matrix_each_model=rbind(posi_matrix_each_model,t(coor_chr_each_model))
      }
      if(model_num==1){
        Dis_each_model=dist(posi_matrix_each_model)
        index_dist_proximal=which(Dis_each_model<=proximity_threshold)
        proximal_pair=sapply(index_dist_proximal,function(x){rowcol(x,nrow(posi_matrix_each_model))})
        proximal_pair=t(rbind(proximal_pair,Dis_each_model[index_dist_proximal]))
      }
      if(model_num!=1){
        dis=apply(proximal_pair,1,function(x){distance_function(posi_matrix_each_model[x[1],],posi_matrix_each_model[x[2],])  })
        proximal_pair=cbind(proximal_pair,dis)
      }
    }
    
    max_dis_all_model=apply(proximal_pair[,3:ncol(proximal_pair)],1,max)
    proximal_pair_consistent=proximal_pair[which(max_dis_all_model<=proximity_threshold),]
    fraction_of_bead_pairs_consistent_within_distance_thresh=rbind(fraction_of_bead_pairs_consistent_within_distance_thresh,
                                                                   c(file,nrow(proximal_pair_consistent)/nrow(proximal_pair)))
    
    
    index_checklist=sapply(1:20,function(x){sum(chr_bead_num_each_file[1:x])})
    index_checklist=c(0,index_checklist)
    
    chr_assign=t(apply(proximal_pair_consistent,1,function(each){
      chr_1=min(which(each[1]<=index_checklist))-1
      chr_2=min(which(each[2]<=index_checklist))-1
      posi_1=each[1]-index_checklist[chr_1]
      posi_2=each[2]-index_checklist[chr_2]
      return(c(chr_1,posi_1,chr_2,posi_2))
    } ))
    
    proximal_pair_consistent=cbind(proximal_pair_consistent,chr_assign)
    proximal_pair_consistent=cbind(proximal_pair_consistent,apply(proximal_pair_consistent[,3:(3+model_num_total-1)],1,mean))
    ncolumn=ncol(proximal_pair_consistent)
    proximal_pair_consistent=proximal_pair_consistent[,c((ncolumn-4):ncolumn,3:(3+model_num_total-1),1:2)]
    ######  pairing in both direction
    proximal_pair_consistent=rbind(proximal_pair_consistent,proximal_pair_consistent[,c(3:4,1:2,5:(ncolumn-2),ncolumn,ncolumn-1)])
    colnames(proximal_pair_consistent)[6]="dis"
    colnames(proximal_pair_consistent)=c("chr_1","bead_label_1","chr_2","bead_label_2","average_distance",colnames(proximal_pair_consistent)[6:(5+model_num_total)],"long_index_1","long_index_2")
    k=k+1
    save(proximal_pair_consistent,file=paste(PATH_output,"/Proximal_pairs_",proximity_threshold,"_",file,".RData",sep=""))
  }
  #print(fraction_of_bead_pairs_consistent_within_distance_thresh)
}

######### below is a general function to map genome positions to bead. 
######## The input should contain the first 3 column of chr; start; end. a bed file R table should be fine.
Get_Bead_from_position=function(input_table,PATH_nuc,PATH_output,reso=1e5,Exact_mapping,chr_format,mode_save=TRUE) {
  library("rhdf5")
  q<-dir(PATH_nuc)
  GenBead_list=vector(mode='list',length=length(q))
  names(GenBead_list)=q
  
  if(chr_format=="chr"){
    chr_num_array=sapply(input_table[,1],function(x){unlist(strsplit(x,split='r'))[2]})
    chr_num_array[which(chr_num_array=='X')]=as.character(20)
    chr_num_absent=which(chr_num_array %in% as.character(1:20)==FALSE)
    chr_num_array[chr_num_absent]=as.character(0)
    chr_num_array=as.numeric(chr_num_array)
    GenTable=cbind(chr_num_array,input_table[,2:3],1:nrow(input_table))
  }
  if(chr_format=="number"){
    GenTable=cbind(input_table[,1:3],1:nrow(input_table))
  }
  if(ncol(input_table)>3){
    GenTable=cbind(GenTable,input_table[,4:ncol(input_table)])
  }
  
  k=1
  for (file in q){
    GenBead=c()
    a=h5read(paste(PATH_nuc,file,sep="/"), "/")
    for (chr_num in 1:20){
      # if(chr_num %in% 1:19 ==TRUE ){ chr_num_old_label=as.character(chr_num)}
      # if(chr_num ==20 ){ chr_num_old_label="X"}
      if(chr_num %in% 1:19 ==TRUE ){ chr_num_old_label=as.character(chr_num)}
      if(chr_num ==20 ){ chr_num_old_label="X"}
      GenTable_chr=GenTable[GenTable[,1]==chr_num,]
      posi1d=a[["structures"]][[1]][["particles"]][[chr_num_old_label]][["positions"]]
      start_1d=posi1d[1]
      end_1d=posi1d[length(posi1d)]
      
      #### 
      if(Exact_mapping==FALSE){
        mean_gene_posi=(GenTable_chr[,2]+GenTable_chr[,3])/2
        beads_num=sapply(mean_gene_posi,function(x){
          (x-start_1d) %/% reso+1
        })
        GenTable_chr=cbind(GenTable_chr,beads_num)
        ################  make sure bead_label is not -1 or -2 etc or bigger than the range of posi1d
        GenTable_chr=GenTable_chr[GenTable_chr[,ncol(GenTable_chr)]>0 & GenTable_chr[,ncol(GenTable_chr)]<=length(posi1d),]
        
      }
      
      if(Exact_mapping==TRUE){
        beads_start=sapply(GenTable_chr[,2],function(x){
          (x-start_1d) %/% reso+1
        })
        beads_end=sapply(GenTable_chr[,3],function(x){
          (x-start_1d) %/% reso+1
        })
        table_new=cbind(GenTable_chr,beads_start,beads_end)#### column number in GenTable_chr can be any number, but the first 4 columns are the same
        ## 
        each_fine_mapping=function(each, posi1d){
          n=length(each)
          if(each[(n-1)]==each[n]){return(each[1:(n-1)])}
          if((each[n]-each[n-1])==1){
            each_start_line=c(each[1],each[2],posi1d[each[n-1]+1]-1,each[4:(n-1)])
            each_end_line=c(each[1],posi1d[each[n]],each[3],each[4:(n-2)],each[n])
            return(rbind(each_start_line,each_end_line))
          }
          if((each[n]-each[n-1])>1){
            each_start_line=c(each[1],each[2],posi1d[each[n-1]+1]-1,each[4:(n-1)])
            each_end_line=c(each[1],posi1d[each[n]],each[3],each[4:(n-2)],each[n])
            each_mid=t(sapply(1:(each[n]-each[n-1]-1),function(y){
              c(each[1],posi1d[each[n-1]+y],posi1d[each[n-1]+y+1]-1,each[4:(n-2)],each[n-1]+y)
            }))
            return (rbind(each_start_line,each_mid,each_end_line))
          }
        }
        ################  make sure bead_label is not -1 or -2 etc or bigger than the range of posi1d
        table_new=table_new[table_new[,(ncol(table_new)-1)]>0 & table_new[,ncol(table_new)]<=length(posi1d),]
        
        table_new2=apply(table_new,1,function(table_new_each_row){t(each_fine_mapping(table_new_each_row,posi1d))})
        GenTable_chr=matrix(unlist(table_new2),ncol=(ncol(GenTable_chr)+1),byrow=TRUE)
      }
      
      GenTable_chr=GenTable_chr[GenTable_chr[,ncol(GenTable_chr)]>0 & GenTable_chr[,ncol(GenTable_chr)]<((end_1d-start_1d)/reso),]
      GenBead=rbind(GenBead,GenTable_chr)
    }
    GenBead=GenBead[,c(1:3,ncol(GenBead),4:(ncol(GenBead)-1))]
    colnames(GenBead)=c("chr","start","end","bead_label","row_num_in_input_table",colnames(input_table)[-c(1,2,3)])
    GenBead_list[[k]]=GenBead
    k=k+1
  }
  if(mode_save=='TRUE'){save(GenBead_list,file=PATH_output)}
  return(GenBead_list)
}
###################

###############################
#### compute intermingling scores for each genome region (here per bead in structure)
Get_intermingling_score_each_file=function(proximal_pair_input,input_GenBead_each_file,reso,thresh_cis_long_range,proximity_thresh){
  proximal_pair_input=proximal_pair_input[proximal_pair_input[,5]<proximity_thresh,]
  intermingling_score_table=c()
  for (chr_num in 1:20){
    proximal_pair_chr=proximal_pair_input[proximal_pair_input[,1]==chr_num,]
    input_GenBead_each_chr=input_GenBead_each_file[input_GenBead_each_file[,1]==chr_num,]
    if(nrow(input_GenBead_each_chr)>1){
      proximal_beads=apply(input_GenBead_each_chr,1,function(x){proximal_pair_chr[x[4]==proximal_pair_chr[,2],1:4]})
    }else{proximal_beads=list()}
    
    get_intermingling_score_each_set=function(proximal_beads_each_set,reso,thresh_cis_long_range){
      if(is.null(nrow(proximal_beads_each_set))){intermingling_score= (c(0,0,0))}
      if(is.null(nrow(proximal_beads_each_set))==FALSE){
        intermingling_score=c(length(which(proximal_beads_each_set[,3]!=proximal_beads_each_set[,1])),
                              length(which(proximal_beads_each_set[,3]==proximal_beads_each_set[,1] & abs(proximal_beads_each_set[,2]-proximal_beads_each_set[,4])>=(thresh_cis_long_range/reso) )),nrow(proximal_beads_each_set))
      }
      return(intermingling_score)
    }
    
    if(length(proximal_beads)>0){
      intermingling_score_table_chr=t(sapply(proximal_beads,function(each){get_intermingling_score_each_set(proximal_beads_each_set=each,reso,thresh_cis_long_range)}))
      if(length(intermingling_score_table_chr)!=0){
        intermingling_score_table=rbind(intermingling_score_table,cbind(intermingling_score_table_chr,input_GenBead_each_chr))
      }
    }
  }
  colnames(intermingling_score_table)[1:3]=c('adj_beads_num_trans','adj_beads_num_cis_long_range','all_adj')
  return(intermingling_score_table)
}

###############

#### wrapper function to literature through different files to compute the intermingling scores for each genome loci within individual cells
#### it returns a list, with each element being a table. The first 3 columns: 'adj_beads_num_trans','adj_beads_num_cis_long_range','all_adj'
## then GenBead_list columns show info related to each row in genome_region_table, namely: chr, start, end, bead_label, input_row_num etc
## calculating the cis long range intermingling score is an addition feature (not used in paper), you can specify a distance above which you want to check in thresh_cis_long_range 
Get_intermingling_score_all_files=function(PATH_proximal,PATH_nuc,genome_region_table,PATH_output,reso=1e5,thresh_cis_long_range=3e6,proximity_thresh){
  q=dir(PATH_proximal)
  intermingling_score_table_all_files=vector(mode="list",length=length(q))
  
  chr_num_old_label_array=vector(mode='integer',length=20)
  chr_num_old_label_array[1:19]=sapply(1:19,function(chr_num){paste('chr',chr_num,sep="")})
  chr_num_old_label_array[20]='chrX'
  chr_num=as.numeric(sapply(genome_region_table[,1],function(x){which(chr_num_old_label_array==x)}))
  
  genome_region_table_new=cbind(chr_num,genome_region_table[,2:3])
  GenBead_list=Get_Bead_from_position(input_table=genome_region_table_new,PATH_nuc,PATH_output=NULL,reso=reso,Exact_mapping=FALSE,chr_format='number',mode_save=FALSE) 
  for(file_num in 1:length(q)){	
    input_GenBead_each_file=GenBead_list[[file_num]]
    load(paste(PATH_proximal,q[file_num],sep=""))
    proximal_pair_consistent=as.data.frame(proximal_pair_consistent,stringsAsFactors = FALSE)
    intermingling_score_table_each_file=Get_intermingling_score_each_file(proximal_pair_input=proximal_pair_consistent,input_GenBead_each_file,reso,thresh_cis_long_range,proximity_thresh)
    intermingling_score_table_all_files[[file_num]]=intermingling_score_table_each_file
  }
  save(intermingling_score_table_all_files,file=PATH_output)
  return(intermingling_score_table_all_files)
}
###########################


####### Please run the pipeline below to get intermingling scores for different cells. 
## Step 1: Use function Dis3D_proximal to find beads pairs in spatial proximity. 
## This can be done with a single command:
Dis3D_proximal(PATH_nuc,PATH_output,proximity_threshold,model_num_total=10)

## PATH_nuc is the path to the folder containing the structural files from all cells 
## PATH_output is the folder where you want to store the results of spatial proximity pairs
## proximity_threshold can normally be set as 2.0 
## model_num_total is the number of models in structural simulation through nuc_dynamics 

###### Step 2 get the intermingling scores for genome regions of interest
## 1. provide a bed file for genome regions you want to calculate the intermingling scores on; this is the 'genome_region_table'
## 2. specifying the proximity file location in PATH_proximal (the output path from Dis3D_proximal)
## 3. reso parameter is the resolution of single cell genome structure 
intermingling_score_all_files_ESC=Get_intermingling_score_all_files(PATH_proximal,PATH_nuc,genome_region_table,PATH_output,reso=1e5,proximity_thresh=2)
