library(dplyr)
library(minpack.lm)
library(mgcv)
library(vegan) #load the R package 'vegan'
library(reshape2)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(broom)
library(tidyr)
library(foreach)
library(pheatmap)
library(ggtext)
source('common_function_bacarena_gapseq.R')


########## Figure S1:check the quantity of reactions that are non-expressed but carry flux
check_ratio_expressed<-function(model,rxns_summary_df){
  fba<-optimizeProb(model)
  index_bio1<-which(grepl("bio1", model@react_id))
  index_firstZF<-which(grepl("zhelunF", model@react_id))[1]
  fluxes<-fba@fluxdist@fluxes@x
  fluxes<-fluxes[-(index_bio1:index_firstZF)] # remove filled in reactions for fluxes and rxn id list same time
  rxn_ID<-model@react_id[-(index_bio1:index_firstZF)] # remove filled in reactions
  active_rxn_ID<-rxn_ID[which(fluxes!=0)]
  active_rxn_ID<-active_rxn_ID[!grepl('_zhelunR',active_rxn_ID)]
  #### check if the reaction is expressed
  active_rxn_ID_trunc<-substr(active_rxn_ID, 1, 8)
  active_rxn_ID_trunc <- active_rxn_ID_trunc[!startsWith(active_rxn_ID_trunc, 'EX_')]
  
  expressed_rxn_all_list<-str_split(paste(rxns_summary_df$modelseed_id,collapse = ' '),' ')[[1]]
  expressed_rxn_all_list_uni<-unique(expressed_rxn_all_list)
  
  active_list<-list()
  inactive_list<-list()
  for (rxn in active_rxn_ID_trunc){
    if (any(grepl(rxn,expressed_rxn_all_list_uni))) active_list<-append(active_list,rxn)
    else inactive_list<-append(inactive_list,rxn)
  }
  
  ratio<-length(active_list)/length(active_rxn_ID_trunc)
  return(list(ratio = ratio, inactive_list = inactive_list))
}

model_names<-c("alt_Areni","bin_16","CCDehalococcoidia_1","CCLitorimonas_sp_1","CCGammaproteobacteria_9","CCNitrospirales_2","CCPseudohongiella_sp_1","bin_6","bin_22","CCAlphaproteobacteria_3","CCAlphaproteobacteria_7","CCGammaproteobacteria_3","CCNitrosopumilus_sp_1","CCRhizobiaceae_merged","CCThermodesulfobacteriota_1")

f_list<-list() # active rxns ratio
active_withoutexpression<-list() # inactive rxns id without expression support

for (model in model_names){
  rxns_summary_df<-read.csv(paste0('summary_TPM',model,'refined.tsv'),sep='\t',quote="",header=TRUE)
  m<-readRDS(paste('CC_symbionts/',model,'_draft_adapt_2024_10_12_expressed_DeadEndFixed_checked.RDS',sep=''))
  
  result<-check_ratio_expressed(m,rxns_summary_df)
  f_list[[model]]<-result[1]
  active_withoutexpression[[model]]<-result[2]
}
percentage_expressed_df_individual<-data.frame(Model_name=names(f_list),Percentage=unlist(f_list),row.names = NULL)
percentage_expressed_df_individual$Model_name<-sapply(percentage_expressed_df_individual$Model_name,rename_MM)
percentage_expressed_df_individual$Model_name<-sapply(percentage_expressed_df_individual$Model_name,function(x) paste0('MM.',x))
cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(percentage_expressed_df_individual$Model_name)))

ggplot(percentage_expressed_df_individual, aes(x = Model_name, y = Percentage * 100, fill = Model_name)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 14) +
  labs(
    y = "Percentage (%)",
    x = "Model name"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),
    legend.position = "none"
  )

########## check the quantity of reactions that are non-expressed but carry flux in the simulatiom
check_ratio_expressed<-function(m,rxns_summary_df){
  index_bio1<-which(grepl("bio1", m@model@react_id))
  index_firstZF<-which(grepl("zhelunF", m@model@react_id))[1]
  fluxes<-m@fbasol$fluxes
  fluxes<-fluxes[-(index_bio1:index_firstZF)] # remove filled in reactions for fluxes and rxn id list same time
  rxn_ID<-m@model@react_id[-(index_bio1:index_firstZF)] # remove filled in reactions
  active_rxn_ID<-rxn_ID[which(fluxes!=0)]
  active_rxn_ID<-active_rxn_ID[!grepl('_zhelunR',active_rxn_ID)]
  #### check if the reaction is expressed
  active_rxn_ID_trunc<-substr(active_rxn_ID, 1, 8)
  active_rxn_ID_trunc <- active_rxn_ID_trunc[!startsWith(active_rxn_ID_trunc, 'EX_')]
  
  expressed_rxn_all_list<-str_split(paste(rxns_summary_df$modelseed_id,collapse = ' '),' ')[[1]]
  expressed_rxn_all_list_uni<-unique(expressed_rxn_all_list)
  
  active_list<-list()
  inactive_list<-list()
  for (rxn in active_rxn_ID_trunc){
    if (any(grepl(rxn,expressed_rxn_all_list_uni))) active_list<-append(active_list,rxn)
    else inactive_list<-append(inactive_list,rxn)
  }
  
  ratio<-length(active_list)/length(active_rxn_ID_trunc)
  return(list(ratio = ratio, inactive_list = inactive_list))
}

#### check the percentage (same as last section) for the community simulation
cycle_range<-110:150
all_results <- list()

for (i in c(1,2,3)){
  replicate_folder<-paste0('CC_symbionts/optimized_',i)
  for (cycle in cycle_range){
    comu<-readRDS(paste0(replicate_folder,'/BacArena_Cymbastella_0121_opt_ex_400grids_mineral_sw_',cycle,'.RDS'))
    f_list<-list() # active rxns ratio
    active_withoutexpression<-list() # inactive rxns id without expression support
    for (MODEL in model_names){
      rxns_summary_df<-read.csv(paste0('summary_TPM',MODEL,'refined.tsv'),sep='\t',quote="",header=TRUE)
      
      m<-comu@specs[[MODEL]]
      result<-check_ratio_expressed(m,rxns_summary_df)
      f_list[[MODEL]]<-result[1]
      active_withoutexpression[[MODEL]]<-result[2]
    }
    print(active_withoutexpression)
    print(unlist(f_list))
    df<-data.frame(Model_name=names(f_list),Percentage=unlist(f_list))
    df$Model_name<-sapply(df$Model_name,rename_MM)
    df$Model_name<-sapply(df$Model_name,function(x) paste0('MM.',x))
    df$replicate <- i
    df$cycle <- cycle
    
    all_results[[length(all_results) + 1]] <- df
  }
}
percentage_expressed_df_simulation <- do.call(rbind, all_results)
avg_percentage_expressed_df_simulation <- percentage_expressed_df_simulation %>%
  group_by(Model_name) %>%
  summarise(
    mean_value = mean(Percentage, na.rm = TRUE),
    sd_value   = sd(Percentage, na.rm = TRUE)
  )
cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(avg_percentage_expressed_df_simulation$Model_name)))

ggplot(avg_percentage_expressed_df_simulation, aes(x = Model_name, y = mean_value * 100, fill = Model_name)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = (mean_value - sd_value) * 100, ymax = (mean_value + sd_value) * 100), 
                width = 0.2, color = "black") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 14) +
  labs(
    y = "Percentage (%)",
    x = "Model name"
  ) +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
    # Tells ggplot to interpret markdown syntax (like *text*) on the Y-axis
    legend.position = "none"
  )

##### combine and make each model two bars ####
percentage_expressed_df_combined <- merge(
  percentage_expressed_df_individual,
  avg_percentage_expressed_df_simulation,
  by = "Model_name"
)
colnames(percentage_expressed_df_combined)<-c('Model','Individual','Simulation','Simulation_SD')

plot_df <- percentage_expressed_df_combined %>%
  pivot_longer(
    cols = c(Individual, Simulation),
    names_to = "Type",
    values_to = "Percentage"
  )
ggplot(plot_df, aes(x = Model, y = Percentage * 100, fill = Type)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    data = subset(plot_df, Type == "Simulation"),
    aes(
      ymin = (Percentage - Simulation_SD) * 100,
      ymax = (Percentage + Simulation_SD) * 100
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Model",
    y = "Percentage (%)",
    fill = ""
  ) +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )
############### PCoA plot ##############
df <- read.csv(paste(getwd,"/CC_symbionts/sum.csv", sep = ''), header = T,sep = ',')
dat       <- sqrt(df[ , 2:16])                     # square-root transform
bray      <- vegdist(dat, method = "bray")
sorensen<-vegdist(dat, method = "bray",binary=TRUE)
dist_vec  <- as.matrix(bray)[1, ]                  # distances to real data (row 1)

# --- build plotting data-frame ---------------------------------------------
plot_df <- tibble(
  medium = df$medium,                              # factor levels
  BC     = dist_vec,
  SR=as.matrix(sorensen)[1, ] 
) %>%
  filter(medium != "real_data")

sr_df<-plot_df %>% group_by(medium)%>% summarise(Average=mean(SR))
plot_df <- plot_df %>%
  mutate(medium = reorder(medium, BC, FUN = median))

# --- draw the boxplot -------------------------------------------------------
p<-ggplot(plot_df, aes(x = medium, y = BC)) +
  # Show the three individual data points
  geom_jitter(
    width = 0.10,
    size = 3,
    alpha = 0.8
  ) +
  
  # Mean
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 4,
    size = 6
  ) +
  
  # Mean +/- SD
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.18,
    linewidth = 0.8
  ) +
  theme_bw() +
  theme(
    axis.title.x  = element_text(size = 15),
    axis.title.y  = element_text(size = 12),
    axis.text.x   = element_blank(),
    axis.text.y   = element_text(size = 12)
  ) +
  labs(x = "Nutrient Condition", y = "Bray–Curtis Dissimilarity")

p

###################################################################
############# crossfeeding plot ############
excluded_cpds<-c('H2O','Biomass','H+','O2','CO2','Hydrogen','Photon','Phosphate','Pyrophosphate','Sodium','Fe2+','Fe3+','Mg','Nickel','K+','Cobalt','Cl-','Ca2+','Cu2+','Zn2+','Mn2+','NAD','NADP','CoA')
object_auxo_cpd<-c('Spermidine','Putrescine','Heme','L-methionine','L-proline','Pantothenic acid','Thiamin')
object_auxo_cpd<-c('Nitrite','Nitric oxide','Nitrate','Nitrous oxide','Ammonium','N2')
object_auxo_cpd<-c('D-Glucosamine','D-Glucose','Glycerol','D-Fructose','Trehalose','Betaine','Succinate','Xylose','D-Mannose','L-Proline','Xylitol','L-Methionine','Choline','Galactose','Maltose','Stachyose')
object_auxo_cpd<-df_major$Metabolites
object_auxo_cpd<-c('Thiamin','Pantothenic acid','Riboflavin','L-Lysine','L-Cysteine','L-Aspartate','L-Serine','L-Threonine','L-Valine','L-Isoleucine','L-Phenylalanine')


exchange_file_1<-read.csv(paste(getwd,"/CC_symbionts/opt_110_150_model_cpd_1.tsv", sep = ''), header = T,sep = '\t')%>%
  filter(met %in% object_auxo_cpd)%>%mutate(prod.flux = prod.flux / 41,cons.flux = cons.flux / 41)
exchange_file_2<-read.csv(paste(getwd,"/CC_symbionts/opt_110_150_model_cpd_2.tsv", sep = ''), header = T,sep = '\t')%>%
  filter(met %in% object_auxo_cpd)%>%mutate(prod.flux = prod.flux / 41,cons.flux = cons.flux / 41)
exchange_file_3<-read.csv(paste(getwd,"/CC_symbionts/opt_110_150_model_cpd_3.tsv", sep = ''), header = T,sep = '\t')%>%
  filter(met %in% object_auxo_cpd)%>%mutate(prod.flux = prod.flux / 41,cons.flux = cons.flux / 41)

exchange_file_1$prod.flux<-abs(exchange_file_1$cons.flux)
exchange_file_2$prod.flux<-abs(exchange_file_2$cons.flux)
exchange_file_3$prod.flux<-abs(exchange_file_3$cons.flux)

corrected_df=make_cf_table(exchange_file_1,exchange_file_2,exchange_file_3)

# corrected_df$prod<-sapply(corrected_df$prod,function(x) rename_MM(x))
# corrected_df$cons<-sapply(corrected_df$cons,function(x) rename_MM(x))
write.table(corrected_df,'cf_CC_110_150.tsv',sep='\t',row.names = FALSE,quote=FALSE)

auxotroph_env<-new.env()
## need this order 
object_auxo_cpd<-c('Spermidine','Putrescine','Heme','L-methionine','L-proline','Pantothenic acid','Thiamin')
auxotroph_env$name<-object_auxo_cpd
auxotroph_env$auxo_models<-list(c('bin_16', 'bin_6', 'CCPseudohongiella_sp_1', 'alt_Areni', 'bin_22', 'CCAlphaproteobacteria_3', 'CCAlphaproteobacteria_7', 'CCDehalococcoidia_1', 'CCGammaproteobacteria_3', 'CCGammaproteobacteria_9', 'CCNitrospirales_2', 'CCRhizobiaceae_merged', 'CCThermodesulfobacteriota_1'),
                                c('bin_16', 'CCRhizobiaceae_merged', 'CCGammaproteobacteria_9', 'CCAlphaproteobacteria_7'),
                                c('CCAlphaproteobacteria_7' ), c('bin_22','CCNitrosopumilus_sp_1'),c('CCNitrosopumilus_sp_1'),
                                c('bin_6', 'bin_22', 'CCLitorimonas_sp_1', 'CCThermodesulfobacteriota_1'),
                                c('CCNitrospirales_2', 'CCRhizobiaceae_merged', 'bin_22', 'CCLitorimonas_sp_1', 'CCAlphaproteobacteria_3', 'CCAlphaproteobacteria_7', 'CCGammaproteobacteria_3', 'CCThermodesulfobacteriota_1'))
# final_df need to be generated by the script section get interaction tables for cytoscape first
final_df$prod<-sapply(final_df$prod,function(x) rename_MM(x))
final_df$cons<-sapply(final_df$cons,function(x) rename_MM(x))

#####################################################################################################
############# bar plot for uptake of organic compounds in basic medium
metabolites_name<-'maltose,fructose,mannose,trehalose,starchyose,glycerol,galactose,glucose,xylose,arabinose,choline,xylitol,succinate,d-glucosamine,betaine,L-proline,L-methionine'
metabolites<-'cpd00179,cpd00082,cpd00138,cpd00794,cpd01133,cpd00100,cpd00108,cpd00027,cpd00154,cpd00224,cpd00098,cpd00306,cpd00036,cpd00276,cpd00540,cpd00129,cpd00060'
commu_name<-'BacArena_Cymbastella_0121_opt_ex_400grids_mineral_sw_'
commu<-paste0('CC_symbionts/Cymbastella_0121_opt_ex_Cymbastella_polysLOW_Co14_sper0.5x_putre_heme_pro_met0.0001_N10x_max1_expressed0.00001DeadEndFixed_',c(1,2,3),'_death_0.05_rmRate_3_NutPer_100_Inoc_0.95_tstep_1h_difspeed_1.7e-09/',commu_name)
cycle<-110:150
linechart_MetList<-function(input_metabolites,cycle,commu_name){
  metstoCheck<-sapply(strsplit(input_metabolites, ",")[[1]],function(x) paste('EX_',x,'_e0',sep=''))
  metList<-list()# the percentage of consumption
  for (i in cycle){
    arena<-readRDS(paste(commu_name,as.character(i),'.RDS',sep=''))
    for (met_vec in metstoCheck){
      after_simu<-sum(arena@medlist[[2]][[met_vec]])/1440000
      initial_simu<-sum(arena@medlist[[1]][[met_vec]])/1440000
      
      metList[[met_vec]]<-c(metList[[met_vec]],1-after_simu/initial_simu)### adjusted to medium concentration (1440000=900x1600), 1 mM in medium equal to 900 in simulation
      #metList[[met_vec]]<-c(metList[[met_vec]],after_simu)### adjusted to medium concentration (1440000=900x1600), 1 mM in medium equal to 900 in simulation
    }
  }
  return(metList)
}

metList_1<-linechart_MetList(metabolites,cycle,commu[1])
df_1<-apply(as.data.frame(metList_1),2,function(x) sum(x)*100 / length(cycle))
metList_2<-linechart_MetList(metabolites,cycle,commu[2])
df_2<-apply(as.data.frame(metList_2),2,function(x) sum(x)*100 / length(cycle))
metList_3<-linechart_MetList(metabolites,cycle,commu[3])
df_3<-apply(as.data.frame(metList_3),2,function(x) sum(x)*100 / length(cycle))

combined_df <- rbind(df_1, df_2, df_3)
result <- data.frame(
  mean = apply(combined_df, 2, mean),
  sd   = apply(combined_df, 2, sd)
)
write.table(result,'met_stats_summary110_150.tsv',sep='\t',row.names = TRUE,quote=FALSE)

data<-read.csv('CC_symbionts/met_stats_summary110_150.tsv',sep='\t',header=TRUE,quote="")
names(data)[1]<-'Compounds'
data$Compounds<-sapply(data$Compounds,switch_metseedid2metname)
data$Compounds<-sub('GLUM','D-Glucosamine',data$Compounds)
scaling_factors <- c(
  "L-Methionine" = 0.1,
  "L-Proline"    = 0.1,
  "Putrescine"   = 0.002,
  "Spermidine"   = 0.0005,
  'Heme'         =0.0005
)

# Loop through each compound and apply the scaling
for (compound in names(scaling_factors)) {
  idx <- which(data$Compounds == compound)
  factor <- scaling_factors[compound]
  data$mean[idx] <- data$mean[idx] * factor
  data$sd[idx]   <- data$sd[idx] * factor
}

data<-data[order(data$mean),]
data$Compounds <- factor(data$Compounds, levels = data$Compounds)  # keep sorted order
data$mean<-data$mean*0.08 # remove the percentage and times 8 to make unit consistent with cf plot
data$sd<-data$sd*0.08
# Draw the bar plot --------------------------------------------------------
ggplot(data, aes(x = mean, y = Compounds)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), height = 0.3) +
  labs(
    x     = "Mean exchange flux (fmol/h)",
    y     = "Compounds"
  ) + theme_minimal(base_size = 12)


######## make the stacked bar for uptake of organic compounds in basic medium ######
#### get corrected_df from the script crossfeeding_plot.R

corrected_df_refined <- corrected_df %>%
  filter(prod == 'Environment'|cons=='Environment') %>%
  select(1:5) %>%
  arrange(met, desc(mean_flux_prod)) %>%
  group_by(met) %>%
  mutate(
    total_flux = sum(mean_flux_prod),
    mean_flux_percent = mean_flux_prod / total_flux * 100,
    cumulative_percent = cumsum(mean_flux_percent),
    sd_prod_percent = sd_prod / total_flux * 100,
    ymin = cumulative_percent - sd_prod_percent,
    ymax = cumulative_percent + sd_prod_percent,
    entry_id = row_number()
  ) %>%
  filter(cumulative_percent <= 80 | entry_id <= 2)
corrected_df_refined$cons<-sapply(corrected_df_refined$cons,rename_MM)

### group_by met and then reorder the df based on alphabetical order to make the error match the stacked bar
corrected_df_refined$cons <- factor(
  corrected_df_refined$cons,
  levels = rev(levels(as.factor(corrected_df_refined$cons)))
)
corrected_df_refined_reorder <- corrected_df_refined %>%
  arrange(met, cons) %>%
  mutate(
    cumulative_percent = rev(cumsum(rev(mean_flux_percent))),
    ymin_percent = cumulative_percent - sd_prod_percent,
    ymax_percent = cumulative_percent + sd_prod_percent)

#### make the percent times absolute flux obtained from the uptake bar plot
corrected_df_refined_reorder <- corrected_df_refined_reorder %>%
  left_join(data, by = c("met" = "Compounds")) %>%
  mutate(
    mean_flux_percent     = mean_flux_percent     * mean,
    sd_prod_percent            = sd_prod_percent       * mean,
    cumulative_percent  = cumulative_percent  * mean,
    ymin               = ymin_percent               * mean,
    ymax               = ymax_percent               * mean
  ) 

write.table(corrected_df_refined_reorder,'met_stats_110_150.tsv',sep='\t',row.names = TRUE,quote=FALSE)

corrected_df_refined_reorder$met<-sapply(corrected_df_refined_reorder$met,correct_capital)

uptake_plot<-ggplot(corrected_df_refined_reorder, aes(x = met, y = mean_flux_percent, fill = cons)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),

    width = 0.3,
    color = "black"
  ) +    scale_fill_manual(values = c("green","#CC00FF","#660033","#CCFF00","blue","pink","#FF9900","yellow","#669933","#33CCFF","#99FFCC","#996600","996666","red",
                                        "#006666","#999999","#CCCCCC","black")) +
  labs(
    x=NULL,
    y = "Mean exchange flux (fmol/h)",
    fill = "Organisms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3),legend.text = element_markdown(size = 8)
  )

# for the release stack bar
ggplot(corrected_df_refined_reorder, aes(x = met, y = mean_flux_percent, fill = prod)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),

    width = 0.3,
    color = "black"
  ) +    scale_fill_manual(values = c("#006666","#FF9900","#CC00FF","#660033","#CCFF00","#CC00FF","#CCFF00","blue","pink","yellow","#669933","#33CCFF","#99FFCC","#996600","996666","red",
                                      "#999999","#CCCCCC","black")) +
  labs(
    x=NULL,
    y = "Mean exchange flux (fmol/h)",
    fill = "Organisms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3)
  )

ggsave("Fig2.tiff", 
       plot = uptake_plot, 
       device = "tiff", 
       width = 200,       # Width in mm (ISME standard full-page width)
       height = 150,      # Adjust height based on your number of metabolites
       units = "mm", 
       dpi = 900,         # This sets your ppi/dpi
       compression = "lzw") # Reduces file size without losing quality

###########################################################################
########### objectively check if the community is stable 
file_1<-read.csv('CC_symbionts/Cymbastella_polys_Co14_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_max1_expressed0.00001DeadEndFixed_1_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')
file_2<-read.csv('CC_symbionts/Cymbastella_polys_Co14_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_max1_expressed0.00001DeadEndFixed_2_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')
file_3<-read.csv('CC_symbionts/Cymbastella_polys_Co14_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_max1_expressed0.00001DeadEndFixed_3_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')

# file_1<-read.csv('CC_symbionts/unconstrain/Cymbastella_polys_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_unconstrain1000to1_1_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')
# file_2<-read.csv('CC_symbionts/unconstrain/Cymbastella_polys_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_unconstrain1000to1_2_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')
# file_3<-read.csv('CC_symbionts/unconstrain/Cymbastella_polys_sper0.5x_putre_heme_pro_met0.0001_N10x_NEWWWW_unconstrain1000to1_3_rm_rate_3_difspeed_1.7e-09_150_0.95.csv',sep=',')

get_dis_r_a<-function(file_1,file_2,file_3){
  ### get the relative abundance of each population in each cycle ####
  relat_abund_df_1<-apply(file_1[,c(-1,-which(colnames(file_1) == 'last_ten_sum_rela'))],2,get_relat_abund)
  relat_abund_df_2<-apply(file_2[,c(-1,-which(colnames(file_2) == 'last_ten_sum_rela'))],2,get_relat_abund)
  relat_abund_df_3<-apply(file_3[,c(-1,-which(colnames(file_3) == 'last_ten_sum_rela'))],2,get_relat_abund)
  
  cycle_count<-ncol(file_1)-2 # get the amount of cycles
  
  relat_abund_df_average<-(relat_abund_df_1+relat_abund_df_2+relat_abund_df_3)/3
  
  nature_relative_abundance<-relat_abund_df_average[,1]
  
  #### switch the relative abundance of the value to the dissimilarities between the simulated value and natural state ####
  ##### use the log to represent the value, if the r_a is equal to nature r_a, the dissimilarities should be 0 (log1).
  dissimilarities_relat_abund_df_ave <- apply(relat_abund_df_average, 2, function(x) {
    # Check for the cell count that are equal to zero and replace them with NA
    log_values <- ifelse(x == 0, NA, log(x/nature_relative_abundance,base = 10))
    return(log_values)
  })
  rownames(dissimilarities_relat_abund_df_ave)<-file_1$X
  colnames(dissimilarities_relat_abund_df_ave)<-0:cycle_count
  max_abs_value <- max(abs(c(min(dissimilarities_relat_abund_df_ave, na.rm = TRUE), 
                             max(dissimilarities_relat_abund_df_ave, na.rm = TRUE))))
  
  # Create a color palette
  color_palette <- colorRampPalette(c("firebrick3", "white", "navyblue"))(100)
  
  # Define breaks, ensuring 0 is exactly in the middle and mapped to white
  breaks <- seq(-max_abs_value, max_abs_value, length.out = length(color_palette))
  results<-list(df=dissimilarities_relat_abund_df_ave,breaks=breaks)
  return(results)
}

df<-apply(file_1[,c(-1,-which(colnames(file_1) == 'last_ten_sum_rela'))],2,get_relat_abund)
df_2<-apply(file_2[,c(-1,-which(colnames(file_2) == 'last_ten_sum_rela'))],2,get_relat_abund)
df_3<-apply(file_3[,c(-1,-which(colnames(file_3) == 'last_ten_sum_rela'))],2,get_relat_abund)
c_mean<-list()
c_sd<-list()
cycle_count<-ncol(file_1)-2
cycle_countfromzero<-ncol(file_1)-2

for (i in 1:cycle_count){
  all_values <- cbind(as.numeric(df[,i]), as.numeric(df_2[,i]), as.numeric(df_3[,i]))
  col_mean <- rowMeans(all_values, na.rm = TRUE)
  col_sd <- apply(all_values, 1, sd, na.rm = TRUE)
  c_mean[[i]]<-col_mean
  c_sd[[i]]<-col_sd
}
cycle<-1:cycle_count 
d_mean<-as.data.frame(c_mean)
d_sd <- as.data.frame(c_sd)

colnames(d_mean) <- colnames(df[,1:cycle_count])
rownames(d_mean) <- file_1$SampleID
colnames(d_sd) <- colnames(df[,1:cycle_count])
rownames(d_sd) <- file_1$SampleID
x_mean <- data.frame(cycle, t(d_mean))
x_sd <- data.frame(cycle, t(d_sd))

data_mean <- reshape2::melt(x_mean, id.vars ='cycle')
data_sd <- reshape2::melt(x_sd, id.vars ='cycle')
# Merge mean and sd data frames by 'cycle' and 'variable'
final_data <- merge(data_mean, data_sd, by = c("cycle", "variable"))
names(final_data) <- c("cycle", "variable", "mean_value", "sd_value")
final_data<-as.data.frame(final_data)

last_cycles <- sort(unique(final_data$cycle), decreasing = TRUE)[1:50]
df_last <- final_data %>% filter(cycle %in% last_cycles)
results <- df_last %>%
  group_by(variable) %>%
  do(tidy(lm(mean_value ~ cycle, data = .))) %>%
  filter(term == "cycle") %>%  # Get only slope terms
  select(variable, slope = estimate, p_value = p.value)
slope_sum<-sum(abs(results$slope))
slope_sum_last<-sum(abs(results$slope))

#final_data<-final_data[which(final_data$cycle>100),]
final_data$variable<-sapply(final_data$variable,rename_MM)

positive_plot<-ggplot(final_data, aes(x = cycle, y = mean_value, color = variable)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2) +
  scale_colour_manual(values = c("green","red","#CCFF00","blue","pink","#FF9900","yellow","#669933","#33CCFF","#99FFCC","#996600","#660033","996666","#CC33FF",
                                 "#006666","#999999","#CCCCCC","#CC00FF","black")) +
  ylab('Relative abundance') + xlab('A')+
  labs(colour = 'MAGs')+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(trans = 'sqrt',expand = c(0, 0))+theme_bw()+theme(panel.grid = element_blank(),legend.text = element_markdown(size = 10))
positive_plot


################### switch multiple medium for pcoa to a matrix ######
# List all CSV files in the folder
file_list <- list.files(path = "Conditions_pcoa", pattern = "\\.csv$", full.names = TRUE)
base_names <- sapply(file_list, function(x) tools::file_path_sans_ext(basename(x)))

# Read all files into a list of data frames
medium_files <- lapply(file_list, function(f) {
  read.csv(f, sep = ",", header = TRUE, quote = "")
})
#excluded_cpds<-c('Al+++','Boric acid','C02','Ca2+','Cl-','Co2+','Cu2+','F-','I-','Fe2+','Fe3+','Mg2+','Mn2+','Nitrate','Nitrite','NH3','Ni2+','O2','Phosphate','K+','H+','sulfide','sulfate','sulfite','H2S','Photon','Thiosulfate','Na+','Zi2+','Zn2+','turanose')
target_cpds<-c('heme','betaine','spermidine','L-methionine','L-proline','putrescine')
filtered_medium<-lapply(medium_files,function(x){
  x<-x[x$Description %in% target_cpds,]
  x[,c('Exchange','Input_mM')]
})

names(filtered_medium)<- base_names
filtered_medium <- lapply(names(filtered_medium), function(nm) {
  df <- filtered_medium[[nm]]
  names(df)[names(df) == "Input_mM"] <- nm
  df
})
merged_df <- Reduce(function(x, y) merge(x, y, by = "Exchange", all = TRUE), filtered_medium)
merged_df[is.na(merged_df)] <- 0
merged_df$Exchange<-sapply(merged_df$Exchange, switch_metseedid2metname)
merged_df<-unique(merged_df)
new_order<-c(levels(p$data$medium))
merged_df <- merged_df[ , c("Exchange", new_order)]

mat <- as.matrix(merged_df[ , -1])          # drop Exchange column
rownames(mat) <- merged_df$Exchange

vals <- sort(unique(c(mat)))## values indicating
col_vec <- c('white', "#FFD92F", "#E5C494", "#A6D854","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
break_vec <- function(x, extend = 1e-8) {
  vals <- sort(unique(as.numeric(x)))        # all distinct levels
  if (length(vals) < 2)
    stop("Need at least two unique values to build breaks.")
  
  mids <- head(vals, -1) + diff(vals) / 2    # mid-points between levels
  c(min(vals) - extend, mids, max(vals) + extend)
}

pheatmap(
  mat,
  color         = col_vec,
  breaks        = break_vec(mat),
  legend_breaks = vals,
  legend_labels = paste(vals, "mM"),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  legend = FALSE
)

df_cf_N<-data.frame()
for (met in levels(as.factor(corrected_df$met))){
  ################# add crossfeeding of this metabolite #############
  df_sub_met<-corrected_df[which(corrected_df$met==met),]
  
  result_met_df<-get_top_interactions(df_sub_met,0.9)
  df_cf_N<-rbind(df_cf_N,result_met_df)
}
df_cf_N$prod <- gsub("\\*", "", df_cf_N$prod)
df_cf_N$cons <- gsub("\\*", "", df_cf_N$cons)

write.table(df_cf_N,'cf_table_nitrogen.tsv',quote=FALSE,sep='\t',row.names=FALSE)


######## major released bar plot (Figure 5) ########
df<-read.csv('CC_symbionts/metabolites_release_110_150_1.tsv',sep='\t',header=TRUE)
df_2<-read.csv('CC_symbionts/metabolites_release_110_150_2.tsv',sep='\t',header=TRUE)
df_3<-read.csv('CC_symbionts/metabolites_release_110_150_3.tsv',sep='\t',header=TRUE)

df_major<-rbind(df,df_2,df_3)%>%group_by(X)%>%summarise(mean_flux=sum(flux)/3,sd=sd(flux))%>%filter(mean_flux>0)
df_major<-df_major %>%mutate(Metabolites = factor(X, levels = X[order(mean_flux)])) # reorder and remove those uninvolved metabolites
df_major <- df_major %>% 
  filter(
    mean_flux > 20,
    !Metabolites %in% c("Sulfate", "Nitrite", "N2", "Phosphate", "Hydrogen peroxide","Ammonium",'Nitric oxide','Hydrogen')
  )

corrected_df_refined <- corrected_df %>%
  filter(cons=='Environment') %>%
  select(1:5) %>%
  arrange(met, desc(mean_flux_prod)) %>%
  group_by(met) %>%
  mutate(
    total_flux = sum(mean_flux_prod),
    mean_flux_percent = mean_flux_prod / total_flux,
    cumulative_percent = cumsum(mean_flux_percent),
    sd_prod_percent = sd_prod / total_flux,
    ymin = cumulative_percent - sd_prod_percent,
    ymax = cumulative_percent + sd_prod_percent,
    entry_id = row_number()
  ) %>%
  filter(
    entry_id == 1                                           # always keep first row
    | cumulative_percent <= 0.8                           # keep rows until we cross 80 %
    | (entry_id == 2 & lag(cumulative_percent) <= 0.8)    # keep 2nd row only if 1st row was still ≤ 0.8
  )
corrected_df_refined$prod<-sapply(corrected_df_refined$prod,rename_MM)

### group_by met and then reorder the df based on alphabetical order to make the error match the stacked bar
corrected_df_refined$prod <- factor(
  corrected_df_refined$prod,
  levels = rev(levels(as.factor(corrected_df_refined$prod)))
)
corrected_df_refined_reorder <- corrected_df_refined %>%
  arrange(met, prod) %>%
  mutate(
    cumulative_percent = rev(cumsum(rev(mean_flux_percent))),
    ymin_percent = cumulative_percent - sd_prod_percent,
    ymax_percent = cumulative_percent + sd_prod_percent)

#### make the percent times absolute flux obtained from the uptake bar plot, divided by 1700 to normalise by biomass
corrected_df_refined_reorder <- corrected_df_refined_reorder %>%
  left_join(df_major, by = c("met" = "Metabolites")) %>%
  mutate(
    mean_flux_nor     = mean_flux_percent     * total_flux/1700,
    sd_prod_nor            = sd_prod_percent       * total_flux/1700,
    ymin               = ymin_percent               * total_flux/1700,
    ymax               = ymax_percent               * total_flux/1700
  ) 

corrected_df_refined_reorder<-corrected_df_refined_reorder[!is.na(corrected_df_refined_reorder$cumulative_percent),]

corrected_df_refined_reorder$met<-sapply(corrected_df_refined_reorder$met,correct_capital)
bar_release_stack<-ggplot(corrected_df_refined_reorder, aes(x = met, y = mean_flux_percent, fill = prod)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(
    aes(
      ymin = pmax(ymin, 0),         # clamp at 0 so sqrt is valid
      ymax = pmax(ymax, 0)
    ),
    
    width = 0.3,
    color = "black"
  ) + scale_fill_manual(values = c("green","#CC00FF","#660033","#CCFF00","blue","pink","#FF9900","yellow","#669933","#33CCFF","#99FFCC","#996600","996666","#CC33FF","red",
                                      "#006666","#999999","#CCCCCC","black")) +
  labs(
    x=NULL,
    y = "Flux (mmol/(g_dw*h))",
    fill = "Organisms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3),legend.text = element_markdown(size = 8)
  )

ggsave("Fig5.tiff", 
       plot = bar_release_stack, 
       device = "tiff", 
       width = 150,       # Width in mm (ISME standard full-page width)
       height = 140,      # Adjust height based on your number of metabolites
       units = "mm", 
       dpi = 900,         # This sets your ppi/dpi
       compression = "lzw") # Reduces file size without losing quality

##### amino acids and cofactor release bar (Figure 6)######
cofactor_release_stack<-ggplot(corrected_df_refined_reorder, aes(x = met, y = mean_flux_nor, fill = prod)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(
    aes(
      ymin = pmax(ymin, 0),         # clamp at 0 so sqrt is valid
      ymax = pmax(ymax, 0)
    ),
    
    width = 0.3,
    color = "black"
  ) + scale_y_continuous(trans = "sqrt") +scale_fill_manual(values = c("green","#CC00FF","#660033","#CCFF00","blue","pink","#FF9900","yellow","#669933","#33CCFF","#99FFCC","#996600","996666","#CC33FF","red",
                                   "#006666","#999999","#CCCCCC","black")) +
  labs(
    x=NULL,
    y = "Flux (mmol/(g_dw*h))",
    fill = "Organisms"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3),legend.text = element_markdown(size = 8)
  )
ggsave("Fig6.tiff", 
       plot = cofactor_release_stack, 
       device = "tiff", 
       width = 150,       # Width in mm (ISME standard full-page width)
       height = 140,      # Adjust height based on your number of metabolites
       units = "mm", 
       dpi = 900,         # This sets your ppi/dpi
       compression = "lzw") # Reduces file size without losing quality


### update the average relative abundance for the table S2
average_RA <- read.table(
  "Relative_abundance_041124.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  select(1:23)

colnames(average_RA)[1] <- "user_genome"
classification_bac <- read.table(
  "gtdbtk.bac120.summary.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  select(c('user_genome', 'classification'))

classification_arc <- read.table(
  "gtdbtk.ar53.summary.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  select(c('user_genome', 'classification'))

quality <- read.table(
  "quality_report.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
) %>%
  select(c('Name', 'Completeness','Contamination'))
colnames(quality)[1] <- "user_genome"

classification_all <- bind_rows(
  classification_bac,
  classification_arc
)

classification_df <- Reduce(
  function(x, y) merge(x, y, by = "user_genome", all.x = TRUE),
  list(classification_all, quality, average_RA)
)

classification_df<-classification_df[order(classification_df$Average,decreasing = TRUE),]
write.table(
  classification_df,file = "classification_df.tsv",
  sep = "\t",quote = FALSE,row.names = FALSE
)

########## Figure 6 plot (crossfeeding interactions) ############
########## load function from the script figure6_network_phyllospongia.R ######
replicate_ids<-c(1,2,3)
file_list<-lapply(replicate_ids, function(x) read.csv(paste0('CC_symbionts/optimized_',x,'/opt_110_150_model_cpd_',x,'.tsv'),sep='\t',header=TRUE,quote=""))
mets<-c('Pantothenic acid','Thiamin')

cf_summary <- make_cf_summary(
  tables = file_list,
  metabolites = mets,
  replicate_ids = replicate_ids,
  missing_as_zero = TRUE
)

auxotroph_env<-tibble(name=c('Pantothenic acid','Thiamin'),
                      auxo_models=list(c('bin_6', 'bin_22', 'CCLitorimonas_sp_1', 'CCThermodesulfobacteriota_1'),
                                        c('CCNitrospirales_2', 'CCRhizobiaceae_merged', 'bin_22', 'CCLitorimonas_sp_1', 'CCAlphaproteobacteria_3', 'CCAlphaproteobacteria_7', 'CCGammaproteobacteria_3', 'CCThermodesulfobacteriota_1'))
)



auxotroph_table <- auxotroph_env %>%
  unnest_longer(auxo_models) %>%
  transmute(
    met = name,
    organism = auxo_models
  ) %>%
  distinct()

# Relative abundances supplied for the Figure 2 taxa.
# These values are already percentages, so they are displayed directly.

model_name <- str_split("CCPseudohongiella_sp_1,CCNitrosopumilus_sp_1,CCNitrospirales_2,bin_16,alt_Areni,CCDehalococcoidia_1,CCThermodesulfobacteriota_1,CCAlphaproteobacteria_3,CCRhizobiaceae_merged,CCLitorimonas_sp_1,CCGammaproteobacteria_9,CCAlphaproteobacteria_7,bin_6,CCGammaproteobacteria_3,bin_22",',')[[1]]

relative_abundance <- as.numeric(str_split(
  '0.319541441,0.2083005,0.104567925,0.070589324,0.056779118,0.035318557,0.032591548,0.028328552,0.052,0.019914428,0.017762807,0.01628389,0.013686997,0.012525176,0.011789105',','
)[[1]])


relative_abundance_lookup <- tibble(
  model_name = model_name,
  relative_abundance = relative_abundance*100
) %>%
  mutate(
    # Flux/auxotrophy tables use the raw cluster-prefixed model IDs.
    organism =  model_name
  ) %>%
  select(organism, relative_abundance)

provider_top80 <- cf_summary$provider_flux %>%
  keep_top_flux(organism_col = "prod", fraction = 0.80)

consumer_top80 <- cf_summary$consumer_flux %>%
  keep_top_flux(organism_col = "cons", fraction = 0.80)

net_cf <- make_net_flux_tables(
  provider_table = provider_top80,
  consumer_table = cf_summary$consumer_flux,
  environment_name = "Environment",
  missing_as_zero = TRUE
)

# Keep raw IDs here so they still match auxotroph_table and node_metadata.
provider_network_table <- net_cf$provider_flux

auxotroph_consumer_rows <- cf_summary$consumer_flux %>%
  inner_join(
    auxotroph_table,
    by = c("met", "cons" = "organism")
  )

consumer_network_table <- bind_rows(
  consumer_top80,
  auxotroph_consumer_rows
) %>%
  distinct(met, cons, .keep_all = TRUE)

node_metadata_template <- make_node_metadata_template(
  provider_table = provider_network_table,
  consumer_table = consumer_network_table,
  figure2_name_fun = rename_MM
)

node_metadata <- node_metadata_template

# Add/override abundance from the values supplied above while retaining any
# metadata abundance for taxa not included in relative_abundance_lookup.
if (!"relative_abundance" %in% names(node_metadata)) {
  node_metadata$relative_abundance <- NA_real_
}

node_metadata <- node_metadata %>%
  mutate(
    organism = as.character(organism),
    metadata_relative_abundance = suppressWarnings(as.numeric(relative_abundance))
  ) %>%
  select(-relative_abundance) %>%
  left_join(relative_abundance_lookup, by = "organism") %>%
  mutate(
    relative_abundance = coalesce(
      relative_abundance,
      metadata_relative_abundance
    )
  ) %>%
  select(-metadata_relative_abundance)

network <- make_three_column_network(
  provider_table = provider_network_table,
  consumer_table = consumer_network_table,
  metabolites = c('Pantothenic acid','Thiamin'),
  auxotrophs = auxotroph_table,
  node_metadata = node_metadata,
  figure2_name_fun = rename_MM,
  metabolite_colours = DEFAULT_METABOLITE_COLOURS,
  min_mean_flux = 0,
  # Supplied abundance values are already percentages (e.g. 13.70 = 13.70%).
  abundance_multiplier = 1,
  abundance_digits = 1
)

network_plot <- plot_three_column_network(
  network,
  show_flux_labels = TRUE
)

network_plot

ggsave(
  filename = "Figure_crossfeeding.pdf",plot = network_plot,width = 14,height = 10,units = "in",bg = "white"
)
# Optional exports with Figure 2 display names. The underlying analysis tables
# above retain raw IDs for reliable joins.
provider_table_for_export <- add_figure2_names_to_flux_table(
  flux_table = net_cf$provider_flux,
  organism_col = "prod",
  node_metadata = node_metadata
)

consumer_table_for_export <- add_figure2_names_to_flux_table(
  flux_table = net_cf$consumer_flux,
  organism_col = "cons",
  node_metadata = node_metadata
)

write.table(
  provider_table_for_export,
  file = "provider_tables8.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  consumer_table_for_export,
  file = "consumer_tables8.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)