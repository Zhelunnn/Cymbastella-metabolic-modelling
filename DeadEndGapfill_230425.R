### gapfill the network based on expressed reactions
library(stringr)
library(optparse)
library(sybil)
source('/srv/scratch/z5245780/scripts_and_shell_command/Bacarena/common_function_bacarena_gapseq.R')
option_list = list(
  optparse::make_option(c("-m", "--model_name"),help="model name "),
  optparse::make_option(c("-u", "--unfeasible_rxns"),type='character',default=NA,help="unfeasible reactions df "),
  optparse::make_option(c("-r", "--reactions_file"),type='character',default=NA,help="the reaction file for this model"),
  optparse::make_option(c("-o", "--output"), type='character',default=NA,help="output the report")
)
opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

model_name<-opt$model_name
un_rxn_file<-opt$unfeasible_rxns
reactions_file<-opt$reactions_file
output<-opt$output

rxn_file<-read.table(reactions_file,sep='\t',quote="",header=TRUE)
un_rxn<-t(read.table(un_rxn_file,sep='\t'))
m<-readRDS(paste0(model_name,'.RDS'))

# rxn_file<-read.table('CCNitrospirales_2-all-Reactions.tbl',sep='\t',quote="",header=TRUE)
# m<-readRDS('CCNitrospirales_2_draft_adapt_2024_10_12_expressed.RDS')
# un_rxn<-t(read.table('unfeasible_reactions.tsv',sep='\t'))
rxns_curate_db<-read.csv('/srv/scratch/z5245780/DB/custom/modelseed_db_github_formatted4curation.tsv',sep='\t',quote="",header = TRUE)

df_list<-lapply(un_rxn,function(x){
  row_index<-which(grepl(x,rxns_curate_db$Rid))[1]
  rxns_curate_db[row_index,]
})
unfeasible_reactions_df<-do.call(rbind,df_list)

# switch unfeasible_reactions_df to S matrix
switch_df2S<-function(unfeasible_reactions_df){
  all_mets <- unique(unlist(strsplit(paste(unfeasible_reactions_df$met_id, collapse=","), ",")))
  S_unfea <- matrix(0, 
                    nrow = length(all_mets), 
                    ncol = nrow(unfeasible_reactions_df), 
                    dimnames = list(all_mets, unfeasible_reactions_df$Rid))
  
  for (i in seq_len(nrow(unfeasible_reactions_df))) {
    mets <- unlist(strsplit(unfeasible_reactions_df$met_id[i], ","))
    coeffs <- as.numeric(unlist(strsplit(unfeasible_reactions_df$met_Scoef[i], ",")))
    S_unfea[mets, unfeasible_reactions_df$Rid[i]] <- coeffs
  }
  lb <- unfeasible_reactions_df$lowbnd
  ub <- unfeasible_reactions_df$uppbnd
  names(lb) <- unfeasible_reactions_df$Rid
  names(ub) <- unfeasible_reactions_df$Rid
  return(list(S=S_unfea,ub=ub,lb=lb))
}
S_unfea<-switch_df2S(unfeasible_reactions_df)

#### get the dead-ends metabolite from model
S_model<-as.matrix(S(m))

lowb <- m@lowbnd
uppb <- m@uppbnd


rownames(S_model)<-m@met_id
colnames(S_model)<-gsub('_c0','_zhelun',m@react_id)

get_dead_end_metabolite<-function(S,ub,lb){
  
  produced <- rep(0, nrow(S))
  consumed <- rep(0, nrow(S))
  
  for (rxn in 1:ncol(S)) {
    # Check if reaction can run forward
    if (ub[rxn] > 0) {
      produced <- produced + (S[, rxn] > 0)
      consumed <- consumed + (S[, rxn] < 0)
    }
    # Check if reaction can run backward
    if (lb[rxn] < 0) {
      produced <- produced + (S[, rxn] < 0)
      consumed <- consumed + (S[, rxn] > 0)
    }
  }
  # Identify dead-ends
  only_produced <- which(produced > 0 & consumed == 0)
  only_consumed <- which(consumed > 0 & produced == 0)
  dead_end_metabs <- unique(c(only_produced, only_consumed))
  dead_end_ids <- rownames(S)[dead_end_metabs]
  dead_end_ids<-dead_end_ids[!grepl('\\[e0\\]',dead_end_ids)]
  return(dead_end_ids)
}

############## to deal with reversible reactions ########################
# Check if a reaction can consume a metabolite
can_consume <- function(stoich_row, lb, ub) {

  (stoich_row < 0 & ub > 0) | (stoich_row > 0 & lb < 0)
}

# Check if a reaction can produce a metabolite
can_produce <- function(stoich_row, lb, ub) {

  (stoich_row > 0 & ub > 0) | (stoich_row < 0 & lb < 0)
}

#met<-'cpd00766[c0]'
rxns_E<-retrieve_unbalanced_rxns()
excluded_rxn<-c('rxn13894_zhelun','rxn15376_zhelun','rxn15377_zhelun')
excluded_rxn <- c(excluded_rxn,rxns_E)
fail_met_duetoDB<-list()
deadend_in_expressedrxns_all<-intersect(get_dead_end_metabolite(S_unfea$S,S_unfea$ub,S_unfea$lb),get_dead_end_metabolite(S_model,m@uppbnd,m@lowbnd))
recursively_eleminateDeadEnd<-function(deadend_in_expressedrxns,added_rxns=NULL,fail_met_list=NULL,model){
  deadend_in_expressedrxns_ori <- deadend_in_expressedrxns  
  if (is.null(added_rxns)) added_rxns <- list()
  if (is.null(fail_met_list)) fail_met_list <- list()
  
  for (met in deadend_in_expressedrxns){
    sub_df<-rxns_curate_db[which(grepl(met,rxns_curate_db$met_id, fixed = TRUE)),]
    sub_df_info<-switch_df2S(sub_df)
    S_sub<-sub_df_info$S ## get all possible path using db
    
    stoich_row <- S_model[which(model@met_id == met), ]
    
    if (any(can_produce(stoich_row, lowb, uppb))){ # if the target compound is the product
      
      stoich_row <- as.vector(S_sub[met, ])  # flatten in case it’s a matrix row
      rxn_lb <- sub_df_info$lb            
      rxn_ub <- sub_df_info$ub

      # Logical vector: TRUE if reaction can produce this metabolite
      can_consume_vec <- can_consume(stoich_row, rxn_lb, rxn_ub)
      # Select only the columns (reactions) where this is TRUE:
      consumed_matrix <- S_sub[, which(can_consume_vec), drop = FALSE]
      
      if (ncol(consumed_matrix)==0){
        fail_met_duetoDB<-append(met,fail_met_duetoDB)
        deadend_in_expressedrxns<-deadend_in_expressedrxns[deadend_in_expressedrxns!=met] # remove fixed met from the dead end vector
        next
      }else{
        cat(paste0(ncol(consumed_matrix),' reactions were identified that consume ',met,'\n'))
      }
      candidate_rxns<-list()
      # loop through mets of each rxn (columns)
      for (col in 1:ncol(consumed_matrix)){
        stoich_vector <- consumed_matrix[, col]
        
        # Select compounds where the stoichiometry sign is opposite to met's sign
        products <- rownames(consumed_matrix)[(stoich_vector * consumed_matrix[met, col] < 0)]
        rxn_name<-colnames(consumed_matrix)[col]
        
        ### check if all product in this reaction can be consumed by the model
        consumed_elsewhere <- sapply(products, function(prod_met) {
          tryCatch({
            any(can_consume(S_model[prod_met, ], model@lowbnd, model@uppbnd))
          },error=function(e){
            FALSE 
          })
        })
        if (all(consumed_elsewhere)) {
          candidate_rxns<-append(rxn_name,candidate_rxns)
          cat(paste0(rxn_name,' could be the candidate consuming reaction for dead end ',met,'\n'))
          deadend_in_expressedrxns<-deadend_in_expressedrxns[deadend_in_expressedrxns!=met] # remove fixed met from the dead end vector
        }
      }
      
    }else if(any(can_consume(stoich_row, lowb, uppb))){
      
      stoich_row <- as.vector(S_sub[met, ])
      rxn_lb <- sub_df_info$lb            
      rxn_ub <- sub_df_info$ub

      can_produce_vec <- can_produce(stoich_row,rxn_lb,rxn_ub)
      produced_matrix <- S_sub[, which(can_produce_vec), drop = FALSE]
      
      candidate_rxns<-list()
      if (ncol(produced_matrix)==0){
        fail_met_duetoDB<-append(met,fail_met_duetoDB)
        deadend_in_expressedrxns<-deadend_in_expressedrxns[deadend_in_expressedrxns!=met] # remove fixed met from the dead end vector
        
        next
      }else{
        cat(paste0(ncol(produced_matrix),' reactions were identified that produce ',met,'\n'))
      }
      
      # loop through mets of each rxn (columns)
      for (col in 1:ncol(produced_matrix)){
        stoich_vector <- produced_matrix[, col]
        
        # Select compounds where the stoichiometry sign is opposite to met's sign
        substrates <- rownames(produced_matrix)[(stoich_vector * produced_matrix[met, col] < 0)]
        rxn_name<-colnames(produced_matrix)[col]
        
        ### check if all product in this reaction can be consumed by the model
        produced_elsewhere <- sapply(substrates, function(sub_met) {
          tryCatch({
            any(can_produce(S_model[sub_met, ], model@lowbnd, model@uppbnd))
          },error=function(e){
            FALSE 
          })
        })
        if (all(produced_elsewhere)) {
          candidate_rxns<-append(rxn_name,candidate_rxns)
          cat(paste0(rxn_name,' could be the candidate producing reaction for dead end ',met,'\n'))
          deadend_in_expressedrxns<-deadend_in_expressedrxns[deadend_in_expressedrxns!=met]# remove fixed met from the dead end vector
        }
      }
    }
    candidate_rxns<-candidate_rxns[!candidate_rxns %in% excluded_rxn]
    if (length(candidate_rxns) == 0){
      cat(paste0(met,' can not be consumed or produced by adding one reaction in the model\n'))
      fail_met_list[[length(fail_met_list)+1]]<-met
    }else if(length(candidate_rxns)==1){ # add the reaction directly, if only one solution found
      curated_row<-sub_df[sub_df$Rid==candidate_rxns[[1]],]
      model<-addReact(model,id=curated_row$Rid,        met = str_split(as.character(curated_row$met_id),',')[[1]],
                      Scoef      = as.integer(str_split(as.character(curated_row$met_Scoef),',')[[1]]),
                      reactName  = as.character(curated_row$react_name),
                      reversible = curated_row$react_rev,
                      metName    = str_split(as.character(curated_row$met_name),',')[[1]],
                      metComp    = as.integer(str_split(as.character(curated_row$met_comp),',')[[1]]),
                      lb=as.numeric(curated_row$lowbnd)/1000,ub=as.numeric(curated_row$uppbnd)/1000)
      cat(paste0(candidate_rxns[[1]],' has been added to model for dead end ',met,'\n'))
      added_rxns<-append(candidate_rxns[[1]],added_rxns)
      
    }else{ # have multiple reactions solution
      # choose the reaction with the highest bitscore
      max_bitscore_list<-list()
      for (i in seq_along(candidate_rxns)) {
        rxn<-gsub('_zhelun','',candidate_rxns[[i]])
        bitscores <- rxn_file[grepl(rxn, rxn_file$dbhit), "bitscore"]
        max_bitscore_list[i] <- if (length(bitscores) > 0) max(bitscores) else -Inf
      }
      # Check for all -Inf (i.e., no valid bitscores)
      if (all(max_bitscore_list == -Inf,na.rm = TRUE)) { ## add the first candidate reaction
        cat("No valid bitscore found for any candidate reaction.First candidate will be added\n")
        selected_indices=1
      }else{
        selected_indices<-which.max(max_bitscore_list)
      }
      for (i in selected_indices){
        curated_row<-sub_df[sub_df$Rid==candidate_rxns[[i]],]
        model<-addReact(model,id=curated_row$Rid,        met = str_split(as.character(curated_row$met_id),',')[[1]],
                        Scoef      = as.integer(str_split(as.character(curated_row$met_Scoef),',')[[1]]),
                        reactName  = as.character(curated_row$react_name),
                        reversible = curated_row$react_rev,
                        metName    = str_split(as.character(curated_row$met_name),',')[[1]],
                        metComp    = as.integer(str_split(as.character(curated_row$met_comp),',')[[1]]),
                        lb=as.numeric(curated_row$lowbnd)/1000,ub=as.numeric(curated_row$uppbnd)/1000)
        cat(paste0(candidate_rxns[[i]],' has been added to model for dead end ',met,'\n'))
        added_rxns<-append(candidate_rxns[[i]],added_rxns)

      }
    }
  }
  
  if (identical(deadend_in_expressedrxns, deadend_in_expressedrxns_ori)){ ### quit the loop until deadend metabolites can never decrease
    cat('Recursive fix completed')
    saveRDS(model,paste0(model_name,'_DeadEndFixed.RDS'))
    return(list(added_rxns=unique(added_rxns),fail_met_list=unique(fail_met_list)))
  }else{
    cat("Proceeding to the next recursive step.\n")
    recursively_eleminateDeadEnd(deadend_in_expressedrxns,added_rxns=added_rxns,fail_met_list = fail_met_list,model)
  }
}

report<-recursively_eleminateDeadEnd(deadend_in_expressedrxns_all,model=m)
writeLines(c(length(report$added_rxns),"\tadded Reactions:", paste(report$added_rxns,collapse = ' '), 
             "\n",length(report$fail_met_list),"\tUnresolved Metabolites:", paste(report$fail_met_list,collapse = ' ')),
           con = output)
