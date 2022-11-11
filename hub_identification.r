## For any queries related to the pipeline below, please contact Xiaoyan Ma (xm227@cam.ac.uk)

## Use the functions below to identify single cell contact hubs from .ncc files generated from Nuc_process
## Note1: this is designed for single haploid cells, which only have only allele for any given loci
## Note2: If wishing to use this for diploid cells, please make sure you disambiguate the diploid 
##        contact map first before putting it into this pipeline to identify hubs. 
## Note3: Please be cautious with ambiguous contacts. It's better only using 
##        reliable unambiguate contacts. 
## Note4: Please run the R function definitions from line 11 to line 380 first, and then 
##        go to line 380 to set input path and run the two functions
##        (i.e. one for cell specific control determination, another one is the
##        main function to identify hubs) from line 384 to line 403
       

## R package dplyr is required, R version >3.4 is preferred
library('dplyr')

#### get cell specific control for one cell (shuffling only long range contacts, while keeping sequence separation the same)
Cell_specific_control=function(input_candid,ncc_file,start_end_all_chr_table,stripe_width,total_run_num){
  stats_all_run=c()
  for(run_num in 1:total_run_num){
    stats_all_chr=rep(0,10)
    for(chr_num in 1:20){
      size_chr=start_end_all_chr_table[start_end_all_chr_table[,1]==chr_num,3]
      input_candid_chr=input_candid[input_candid[,1]==chr_num,]
      rand_shift=sample(1:size_chr,nrow(input_candid_chr))
      rand_contact=cbind(input_candid_chr[,2]+rand_shift,input_candid_chr[,5]+rand_shift)
      to_substract=unique(c(which(rand_contact[,1]>size_chr),which(rand_contact[,2]>size_chr)))
      rand_contact[to_substract,1]=rand_contact[to_substract,1]-size_chr
      rand_contact[to_substract,2]=rand_contact[to_substract,2]-size_chr
      rand_contact=rand_contact[rand_contact[,1]>0 & rand_contact[,2]>0,]
      rand_contact=rand_contact[sample(1:nrow(rand_contact),nrow(input_candid_chr)/2),] ## this is because we are going to do dup on it later... this is already on a duplicated set. 
      rand_contact=rbind(rand_contact,cbind(rand_contact[,2],rand_contact[,1]))
      
      chr_start=start_end_all_chr_table[start_end_all_chr_table[,1]==chr_num,2]
      chr_end=start_end_all_chr_table[start_end_all_chr_table[,1]==chr_num,3]
      se1=seq(chr_start,chr_end-stripe_width*2,stripe_width*2)
      se2=seq(chr_start+stripe_width*2,chr_end,stripe_width*2)
      chr_bin_table=cbind(rep(chr_num,length(se1)),se1,se2)
      
      #### 
      map_to_bin=Get_block_from_position_cell_coords_the_same(input_table=cbind(rep(chr_num,nrow(rand_contact)),rand_contact[,1],rand_contact[,1]),block_table=chr_bin_table[,1:3],chr_format='number',mode_save=FALSE)
      reads_abundence_each_bin=sapply(1:nrow(chr_bin_table),function(x){length(which(map_to_bin[,4]==x))})
      stats_chr=sapply(0:29,function(n){length(which(reads_abundence_each_bin==n))})
      stats_all_chr=stats_all_chr+stats_chr
    }
    stats_all_run=rbind(stats_all_run,stats_all_chr/sum(stats_all_chr))
  }
  return(stats_all_run)
}

### The function below is a more general function, which is used to bin loci along the genome coordinate
### block_table/input table: chr (number), start, end; 
### only chr1-chr19 and chrX are considered. This is for mouse genome. 
## If using it for human genome or some other organism, please modify line 7-10 below
Get_block_from_position_cell_coords_the_same=function(input_table,block_table,PATH_output,chr_format,mode_save,return_all_mapped=FALSE) {
  
  block_table=cbind(block_table,1:nrow(block_table))
  
  GenBead_table=c()
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
  
  GenBead=c()
  
  if(return_all_mapped==FALSE){
    for (chr_num in 1:20){
      block_table_chr=block_table[block_table[,1]==chr_num,]
      GenTable_chr=GenTable[GenTable[,1]==chr_num,] 
      if(mode(block_table_chr)=='list' & length(block_table_chr)>0){block_table_chr=as.data.frame(block_table_chr)}
      if(mode(GenTable_chr)=='list' & length(GenTable_chr)>0){GenTable_chr=as.data.frame(GenTable_chr)}
      if((is.matrix(block_table_chr)==TRUE | is.data.frame(block_table_chr)==TRUE) 
         & (is.matrix(GenTable_chr)==TRUE | is.data.frame(GenTable_chr)==TRUE)){if(nrow(block_table_chr)>=1 &nrow(GenTable_chr)>=1){
           mean_gene_posi=(GenTable_chr[,2]+GenTable_chr[,3])/2
           block_num=unlist(sapply(mean_gene_posi,function(x){
             block_num_each=block_table_chr[which(block_table_chr[,2]<=x & block_table_chr[,3]>x),ncol(block_table)]
             if(length(block_num_each)>0){return(block_num_each[1])}else{return(-1)}
           }))
           GenTable_chr=cbind(GenTable_chr,block_num)
           
           GenBead=rbind(GenBead,GenTable_chr)
         }
      }
    }
  }
  
  if(return_all_mapped==TRUE){
    for (chr_num in 1:20){
      block_table_chr=block_table[block_table[,1]==chr_num,]
      GenTable_chr=GenTable[GenTable[,1]==chr_num,]
      if((is.matrix(block_table_chr)==TRUE | is.data.frame(block_table_chr)==TRUE )
         & (is.matrix(GenTable_chr)==TRUE | is.data.frame(GenTable_chr)==TRUE)){if(nrow(block_table_chr)>1 &nrow(GenTable_chr)>1){
           
           mean_gene_posi=(GenTable_chr[,2]+GenTable_chr[,3])/2
           GenTable_chr_with_block_num=c()
           for(i in 1:length(mean_gene_posi)){
             block_num_each=block_table_chr[which(block_table_chr[,2]<=mean_gene_posi[i] & block_table_chr[,3]> mean_gene_posi[i]),ncol(block_table)]
             if(length(block_num_each)>0){
               each=cbind(matrix(rep(GenTable_chr[i,],length(block_num_each)),ncol=ncol(GenTable_chr),byrow=TRUE),block_num_each)
               GenTable_chr_with_block_num=rbind(GenTable_chr_with_block_num,each)
             }
             if(length(block_num_each)==0){
               each=c(GenTable_chr[i,],-1)
               GenTable_chr_with_block_num=rbind(GenTable_chr_with_block_num,each)
             }
           }
           GenBead=rbind(GenBead,GenTable_chr_with_block_num)
         }
      }
    }
    GenBead=apply(GenBead,2,function(x){as.numeric(x)})
  }
  
  if(length(which(GenTable[,1]==0))>=1){
    GenTable_chr0=GenTable[GenTable[,1]==0,]
    block_num=rep(-1,nrow(GenTable_chr0))
    to_put_in_results=cbind(GenTable_chr0,block_num)
    GenBead=rbind(GenBead,to_put_in_results)
  }
  if(length(GenBead)>0){
    GenBead=GenBead[,c(1:3,ncol(GenBead),4:(ncol(GenBead)-1))]
    colnames(GenBead)=c("chr","start","end","block_label","row_num_in_input_table",colnames(input_table)[-c(1,2,3)])
  }else(GenBead=c())
  GenBead_table=GenBead
  if(mode_save=='TRUE'){save(GenBead_table,file=PATH_output)}
  return(GenBead_table)
}
############



## wrapper function to get cell specific control for all cells (all files in one folder) and to produce an imperical p-value look up table 
####  value: p_look_up_table_all_files: columns are ordered as (file_num, p_0,p_1,p_2,etc.. to p_29); (p_i is the for the number of i long range contacts)
Cell_specific_control_all_files=function(FOLDER_PATH_ncc,start_end_all_chr_table,stripe_width=1e4,total_run_num=100,long_range_threshold=5e4,max_cutoff=5e6){
  q<-dir(FOLDER_PATH_ncc)
  p_look_up_table_all_files=c()
  for(file_num in 1:length(q)){
    file=q[file_num]
    ncc_file=read.table(paste(FOLDER_PATH_ncc,file,sep=''),skip=0,stringsAsFactors=FALSE,header=FALSE)
    
    chr_num_old_label_array=vector(mode='integer',length=20)
    chr_num_old_label_array[1:19]=sapply(1:19,function(chr_num){paste('chr',chr_num,sep="")})
    chr_num_old_label_array[20]='chrX'
    chr_first=as.numeric(sapply(ncc_file[,1],function(x){which(chr_num_old_label_array==x)}))
    chr_second=as.numeric(sapply(ncc_file[,7],function(x){which(chr_num_old_label_array==x)}))
    ncc=cbind(chr_first,ncc_file[,2:3],chr_second,ncc_file[,c(8:9,13)])
    ncc=cbind(ncc,1:nrow(ncc))
    ncc=ncc[is.na(ncc[,1])==FALSE & is.na(ncc[,4])==FALSE,]
    colnames(ncc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
    
    ncc_another_direc=cbind(ncc[,4:6],ncc[,1:3],ncc[,7:8])
    colnames(ncc_another_direc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
    ncc_dup=rbind(ncc,ncc_another_direc)
    ncc_dup=ncc_dup[order(ncc_dup[,1],ncc_dup[,2]),]### it matters ! 
    ncc_dup=ncc_dup[-which(duplicated(ncc_dup[,1:2])==TRUE),] ### just in case 
    
    ncc_dup_cis_long_range=ncc_dup[ncc_dup[,1]==ncc_dup[,4] & abs(ncc_dup[,5]-ncc_dup[,2])>long_range_threshold & abs(ncc_dup[,5]-ncc_dup[,2])<max_cutoff,]
    input_candid=ncc_dup_cis_long_range
    
    stats_control=Cell_specific_control(input_candid,ncc_file,start_end_all_chr_table,stripe_width,total_run_num)
    stats_control_average=apply(stats_control,2,function(x){mean(x)})
    p_look_up_table=c(1,sapply(1:(length(stats_control_average)-1),function(x){1-sum(stats_control_average[1:x])}))
    p_look_up_table_all_files=rbind(p_look_up_table_all_files,c(file_num,p_look_up_table))
  }
  return(p_look_up_table_all_files)
  
}

################

### hub is previously called 'stripe' as an analogy for stripes in population Hi-C contact map, 
### that's why hub is referred to as stripe in the code below.
### side note: chr column in input_candid is integer (chr1-chr19 is 1-19, chrX is labeled as 20); input_candid:3 columns minimum; 
Get_stripe_info=function(input_candid,ncc_file,cell_specific_control,stripe_width=1e4,min_span=1e6,min_num_distinct_points=4,anchor_range_thresh=5e6,merge=TRUE){
  library('dplyr')
  short_range_exc_thresh=5e4
  distinct_point_thresh=short_range_exc_thresh ##
  min_points_num_near_anchor=min_num_distinct_points
  
  chr_num_old_label_array=vector(mode='integer',length=20)
  chr_num_old_label_array[1:19]=sapply(1:19,function(chr_num){paste('chr',chr_num,sep="")})
  chr_num_old_label_array[20]='chrX'
  chr_first=as.numeric(sapply(ncc_file[,1],function(x){which(chr_num_old_label_array==x)}))
  chr_second=as.numeric(sapply(ncc_file[,7],function(x){which(chr_num_old_label_array==x)}))
  ncc=cbind(chr_first,ncc_file[,2:3],chr_second,ncc_file[,c(8:9,13)])
  ncc=cbind(ncc,1:nrow(ncc))
  ncc=ncc[is.na(ncc[,1])==FALSE & is.na(ncc[,4])==FALSE,]
  colnames(ncc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
  
  ncc_another_direc=cbind(ncc[,4:6],ncc[,1:3],ncc[,7:8])
  colnames(ncc_another_direc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
  
  ncc_dup=rbind(ncc,ncc_another_direc)
  ncc_dup=ncc_dup[order(ncc_dup[,1],ncc_dup[,2]),]
  ncc_dup=ncc_dup[-which(duplicated(ncc_dup[,1:2])==TRUE),] 
  
  ncc_dup_cis_short_range_exc=ncc_dup[ncc_dup[,1]==ncc_dup[,4] & abs(ncc_dup[,5]-ncc_dup[,2])>short_range_exc_thresh,]
  ncc_dup_trans=ncc_dup[ncc_dup[,1]!=ncc_dup[,4],]
  
  stripe_anchors=c()
  anchor_contacts_info=c()
  stripe_anchor_num=1
  for(chr_num in 1:length(chr_num_old_label_array)){
    ncc_dup_cis_short_range_exc_chr=ncc_dup_cis_short_range_exc[ncc_dup_cis_short_range_exc[,1]==chr_num,]
    candid_chr=input_candid[input_candid[,1]==chr_num,]  ## it can have >3 columns
    if(nrow(candid_chr)>0){
      for(i in 1:nrow(candid_chr)){
        each_point=as.numeric(candid_chr[i,])
        potential_contacts=ncc_dup_cis_short_range_exc_chr[which(abs(ncc_dup_cis_short_range_exc_chr[,2]-each_point[2])<stripe_width),]
        potential_contacts=potential_contacts[order(potential_contacts[,5]),]
        span=potential_contacts[nrow(potential_contacts),5]-potential_contacts[1,5]
        
        if(span>min_span & nrow(potential_contacts)>=min_num_distinct_points){
          distinct_p_list=potential_contacts[1,5]
          for(k in potential_contacts[2:nrow(potential_contacts),5]){
            if((k-distinct_p_list[length(distinct_p_list)])>distinct_point_thresh){
              distinct_p_list=c(distinct_p_list,k)}
          }
          if(length(distinct_p_list)>=min_num_distinct_points){
            stat=c(nrow(potential_contacts),length(distinct_p_list),span)
            if(length(which(abs(distinct_p_list-each_point[2])<anchor_range_thresh))>=min_points_num_near_anchor){
              is_stripe=1
            }else if (length(which(abs(distinct_p_list-each_point[2])<anchor_range_thresh))==0){
              is_stripe=-1
            }else{is_stripe=0}
            if(length(each_point)>3){
              stripe_anchor_each=c(each_point[1:3],stripe_anchor_num,stat,is_stripe,each_point[4:length(each_point)])
            }
            if(length(each_point)==3){
              stripe_anchor_each=c(each_point[1:3],stripe_anchor_num,stat,is_stripe)
            }
            stripe_anchors=rbind(stripe_anchors,stripe_anchor_each)
            
            anchor_contacts_info=rbind(anchor_contacts_info,
                                       cbind(matrix(rep(stripe_anchor_each[1:8],nrow(potential_contacts)),ncol=8,byrow=TRUE),potential_contacts))
            stripe_anchor_num=stripe_anchor_num+1
          }
        }
      }
    }
  }
  stripe_anchors=as.data.frame(stripe_anchors)
  colnames(stripe_anchors)[1:8]=c('chr','start','end','anchor_num','all_contacts_num','unique_contacts_num','span','is_stripe')
  colnames(anchor_contacts_info)[1:8]=c('chr','start','end','anchor_num','all_contacts_num','unique_contacts_num','span','is_stripe')
  ### note: stripe_anchor_num is not ordered according to genome position
  
  stripe_anchors->stripe_anchors_old
  stripe_anchors=stripe_anchors_old[stripe_anchors_old[,8]==1,]
  stripe_anchors_other_type=stripe_anchors_old[stripe_anchors_old[,8]!=1,]
  stripe_anchors=stripe_anchors[order(stripe_anchors[,1],stripe_anchors[,2]),]
  
  anchor_contacts_info->anchor_contacts_info_old
  anchor_contacts_info=anchor_contacts_info_old[anchor_contacts_info_old[,8]==1,]
  anchor_contacts_info_other_type=anchor_contacts_info_old[anchor_contacts_info_old[,8]!=1,]
  anchor_contacts_info=anchor_contacts_info[order(anchor_contacts_info[,1],anchor_contacts_info[,2],anchor_contacts_info[,9],anchor_contacts_info[,10]),]
  
  ## resolving overlapping stripes
  confirmed_anchors=c()
  for (chr_num in 1:length(chr_num_old_label_array)){
    stripe_anchors_chr=stripe_anchors[stripe_anchors[,1]==chr_num,]
    anchor_contacts_info_chr=anchor_contacts_info[anchor_contacts_info[,1]==chr_num,]
    confirmed_anchor_chr=c()
    covered_posi=0
    i=1
    while(i<=nrow(stripe_anchors_chr)){
      if(stripe_anchors_chr[i,2]<=covered_posi){i=i+1}else{
        anchor_point=stripe_anchors_chr[i,2]
        points=anchor_contacts_info_chr[anchor_contacts_info_chr[,2]==anchor_point,10]
        
        points=points[points>covered_posi]
        
        potential_points_all=points
        for(point_each in points){
          potential_points=anchor_contacts_info_chr[anchor_contacts_info_chr[,2]==point_each,10] 
          if(length(which(potential_points %in% points))>=2){potential_points_all=c(potential_points_all,potential_points)}
        }
        potential_points_all=unique(potential_points_all)
        potential_anchors=potential_points_all[potential_points_all%in% stripe_anchors_chr[,2]]
        potential_anchors=potential_anchors[order(potential_anchors)]
        
        potential_anchors=potential_anchors[potential_anchors>covered_posi]
        
        unique_contact_num=sapply(potential_anchors,function(x){stripe_anchors_chr[stripe_anchors_chr[,2]==x,6]})
        potential_chosen_index=which(unique_contact_num==max(unique_contact_num))
        chosen_index=potential_chosen_index[length(potential_chosen_index) %/%2+1 ]
        chosen_anchor=potential_anchors[chosen_index]
        confirmed_anchor_chr=c(confirmed_anchor_chr,chosen_anchor)
        
        points_within_stripe=anchor_contacts_info_chr[anchor_contacts_info_chr[,2]==chosen_anchor,10]
        if(length(points_within_stripe)>0){  #### just to make sure, normally no exception
          covered_posi=points_within_stripe[length(points_within_stripe)]
        }
        i=i+1
      }
    }
    confirmed_anchors=rbind(confirmed_anchors,cbind(rep(chr_num,length(confirmed_anchor_chr)),confirmed_anchor_chr))
  }
  
  colnames(confirmed_anchors)=c('chr','start')
  confirmed_anchors=as.data.frame(confirmed_anchors)
  anchor_contacts_info=as.data.frame(anchor_contacts_info)
  stripe_anchors=as.data.frame(stripe_anchors)
  
  confirmed_anchor_table=inner_join(confirmed_anchors,stripe_anchors,by=c('chr'='chr','start'='start'))
  confirmed_anchor_table=confirmed_anchor_table[order(confirmed_anchor_table[,1],confirmed_anchor_table[,2]),]
  
  confirmed_anchor_contacts_info=inner_join(confirmed_anchors,anchor_contacts_info,by=c('chr'='chr','start'='start'))
  confirmed_anchor_contacts_info=confirmed_anchor_contacts_info[order(confirmed_anchor_contacts_info[,1],confirmed_anchor_contacts_info[,2],confirmed_anchor_contacts_info[,9],confirmed_anchor_contacts_info[,10]),]
  
  ## further filter stripes according to FDR q values
  #### FDR control (0.05)
  p_value=cbind(confirmed_anchor_contacts_info[,4],sapply(confirmed_anchor_contacts_info[,5],function(contact_num){cell_specific_control[contact_num+1]}))
  p_value=p_value[!duplicated(p_value[,1]),]
  total_bin_num=sum((start_end_all_chr_table[,3]-start_end_all_chr_table[,2])/2/stripe_width)
  p_adjusted=p.adjust(p_value[,2],n=total_bin_num,method='BH')
  anchor_num_chosen=p_value[which(p_adjusted<0.05),1]
  confirmed_anchor_contacts_info=confirmed_anchor_contacts_info[confirmed_anchor_contacts_info[,4] %in% anchor_num_chosen,]
  print(c(file_num,min(confirmed_anchor_contacts_info[,5])))
  
  if(merge){
    return(confirmed_anchor_contacts_info)
  }
  
  if(merge==FALSE){
    return(stripe_anchors_old)
  }
}
####################

## the main function (wrapper function) to identify hubs from all ncc files
## output is 'PATH/output_file_name' where hub contact information is saved
identify_hubs = function(FOLDER_PATH_ncc,cell_specific_control,output,stripe_width=2e4){
  
  long_range_threshold=5e4 ## we only consider long range contacts for hub identification
  
  q<-dir(FOLDER_PATH_ncc)
  stripe_info_list=vector(mode='list',length=length(q))
  
  for(file_num in 1:length(q)){
    file=q[file_num]
    ncc_file=read.table(paste(FOLDER_PATH_ncc,file,sep=''),skip=0,stringsAsFactors=FALSE,header=FALSE)
    
    chr_num_old_label_array=vector(mode='integer',length=20)
    chr_num_old_label_array[1:19]=sapply(1:19,function(chr_num){paste('chr',chr_num,sep="")})
    chr_num_old_label_array[20]='chrX'
    chr_first=as.numeric(sapply(ncc_file[,1],function(x){which(chr_num_old_label_array==x)}))
    chr_second=as.numeric(sapply(ncc_file[,7],function(x){which(chr_num_old_label_array==x)}))
    ncc=cbind(chr_first,ncc_file[,2:3],chr_second,ncc_file[,c(8:9,13)])
    ncc=cbind(ncc,1:nrow(ncc))
    ncc=ncc[is.na(ncc[,1])==FALSE & is.na(ncc[,4])==FALSE,]
    colnames(ncc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
    
    ncc_another_direc=cbind(ncc[,4:6],ncc[,1:3],ncc[,7:8])
    colnames(ncc_another_direc)=c('chr_first','start_first','end_first','chr_second','start_second','end_second','unique_identifier_of_contact','input_ncc_file_row_num')
    ncc_dup=rbind(ncc,ncc_another_direc)
    ncc_dup=ncc_dup[order(ncc_dup[,1],ncc_dup[,2]),]
    ncc_dup=ncc_dup[-which(duplicated(ncc_dup[,1:2])==TRUE),] 
    
    ncc_dup_cis_long_range=ncc_dup[ncc_dup[,1]==ncc_dup[,4] & abs(ncc_dup[,5]-ncc_dup[,2])>long_range_threshold,]
    input_candid=ncc_dup_cis_long_range
    
    cell_specific_control=cell_specific_control_all_files[cell_specific_control_all_files[,1]==file_num,-1]
    stripe_info=Get_stripe_info(input_candid,ncc_file,cell_specific_control,stripe_width=2e4,min_span=1e5,min_num_distinct_points=3,anchor_range_thresh=5e6)
    print(c(nrow(stripe_info),length(unique(stripe_info[,4]))))
    stripe_info_list[[file_num]]=stripe_info
  }
  
  saveRDS(stripe_info_list,file=output,version=2)
}


### below is the pipeline to identify single cell contact hubs
### hub is previously called 'stripe' as an analogy for stripes in population Hi-C contact map. That's why hub is referred to as stripe in the code below.

## below specifies the path to ncc files(path to the folder that contains all ncc files for all cells, which you want to find hubs from)
FOLDER_PATH_ncc='...' 

## start_end_all_chr_table is a table specifying the start and the end of each chromosome, this could be different for different versions of the genome
## for mouse genome, if you want to use the example from github site, please download
##    single_cell_structure_model_chr_start_and_end.csv and save it in the same folder together with other ncc files
## the table provided below is just for example, you can make such a table for your own genome version; 
## it's similar to bed format, except the first column, where chromosome is named 1-19 rather than chr1-chr19. chrX here is labeled as 20. chrY and chrM are not considered
start_end_all_chr_table=read.table(paste(FOLDER_PATH_ncc,'single_cell_structure_model_chr_start_and_end.csv',sep=''),sep=',',stringsAsFactors = FALSE,header = F)

## get cell specific controls within each cell
cell_specific_control_all_files=Cell_specific_control_all_files(FOLDER_PATH_ncc,start_end_all_chr_table,stripe_width=2e4,total_run_num=400)
## if you want to save running time, reduce total_run_num

## run the function below to get all hub contacts from all ncc files
## output is 'PATH/output_file_name' where hub contact information is saved
identify_hubs(FOLDER_PATH_ncc,cell_specific_control,output='...',stripe_width=2e4)
  