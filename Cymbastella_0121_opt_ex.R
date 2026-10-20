library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)
library(parallel)
library(glpkAPI)

option_list = list(
  optparse::make_option(c("-o", "--otu_name"),                 type="character",      default='Cymbastella_0121_opt_ex',       help="OTU_name"),
  optparse::make_option(c("-D", "--infile_BacarenaRDS_diet"),  type="character",
                        help="Diet RDS file for bacarena with full path"),
  optparse::make_option(c("-k", "--keywd"),                    type="character",      help="(Nutirent/date) keyword used in the folder name of output files, defalt: sw_v20210716"),
  optparse::make_option(c("-r", "--rm_rate"),   type="double",       default=0,        help="removeM value, defalt: 0"),
  optparse::make_option(c("-p", "--tstep"),     type="double",       default=1,         help="Time step, defalt: 1h per iteration"),
  
  optparse::make_option(c("-a", "--arena_mn"),              type="double",       default=20,        help="integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),               type="double",       default=0,         help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-i", "--inocc_no"),              type="double",       default=1,         help="inocculum, defalt: 1"),
  optparse::make_option(c("-c", "--cl_no"),                 type="double",       default=1,         help="Number of replicates, defalt: 1"),
  optparse::make_option(c("-s", "--setAllExInf_value"),     type="character",    default='FALSE',   help="setAllExInf, defalt: FALSE"),
  optparse::make_option(c("-g", "--NutConPercent"),         type="double",       default=100,       help="percentage of nutrients, defalt: 100 (unit: %)"),
  optparse::make_option(c("-n", "--auto_num"),              type="double",       default=1,         help="surfix in the simulation RDS file name."),
  optparse::make_option(c("-S", "--difspeed"),              type="double",       default=1.7e-09,   help="A number indicating the diffusion speed (given by number of cells per iteration), defalt: 1.7e-09 cm2 h-1"),
  optparse::make_option(c("-t", "--iter"),                  type="double",       default=5,         help="Number of iteration, defalt:5"));

opt_parser = optparse::OptionParser(option_list=option_list, add_help_option=FALSE);
opt = optparse::parse_args(opt_parser);

otu_name        = opt$otu_name
infile_diet     = opt$infile_BacarenaRDS_diet
keywd           = opt$keywd
grid_no   = opt$arena_mn
death_r   = opt$death_r
Inocc_no  = opt$inocc_no
Cl_no     = opt$cl_no
setAllExInf_value = opt$setAllExInf_value
NutConPercent   = opt$NutConPercent
iter_no   = opt$iter
rm_rate   = opt$rm_rate##########removing may take place between last iter and this iter,check number of cell by simlist[['type]]==n
time_step =  opt$tstep
auto_num  = opt$auto_num

difspeed = opt$difspeed


#####################################################################################################################################################################

############ Katana ############ 
diet <- readRDS(infile_diet)#/srv/scratch/z5245780/SpongeMAGs_Robbins_ISME/all_sponge_DB_bins_over50pctComp/CAR/genome/Bacarena_community/Bacarena_medium/commounity-medium_seawater_0713.RDS
#####################################################################################################################################################################
getwd <- getwd()
getwd
# create a new folder:

new_folder = paste(otu_name,'_',keywd,'_death_',death_r,'_rmRate_',rm_rate,'_NutPer_',NutConPercent,'_Inoc_',Inocc_no, '_tstep_',time_step,'h_difspeed_',difspeed,'/',sep = '')
dir.create(paste(getwd,'/',new_folder,sep = ''))
dir_gs_models='/srv/scratch/z5245780/Cymbastela_concentrica/sanger_sequencing/CC_sanger_final_bins/bacarena/models'


model1 <- readRDS(paste(dir_gs_models,'/alt_Areni_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model2 <- readRDS(paste(dir_gs_models,'/bin_6_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model3 <- readRDS(paste(dir_gs_models,'/bin_16_ctg8312_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model4 <- readRDS(paste(dir_gs_models,'/bin_22_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model5 <- readRDS(paste(dir_gs_models,'/CCAlphaproteobacteria_3_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model6 <- readRDS(paste(dir_gs_models,'/CCAlphaproteobacteria_7_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model7 <- readRDS(paste(dir_gs_models,'/CCDehalococcoidia_1_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model8 <- readRDS(paste(dir_gs_models,'/CCGammaproteobacteria_3_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model9 <- readRDS(paste(dir_gs_models,'/CCGammaproteobacteria_9_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model10 <- readRDS(paste(dir_gs_models,'/CCLitorimonas_sp_1_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model11 <- readRDS(paste(dir_gs_models,'/CCNitrosopumilus_sp_1_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model12 <- readRDS(paste(dir_gs_models,'/CCNitrospirales_2_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model13 <- readRDS(paste(dir_gs_models,'/CCPseudohongiella_sp_1_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model14 <- readRDS(paste(dir_gs_models,'/CCRhizobiaceae_merged_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))
model15 <- readRDS(paste(dir_gs_models,'/CCThermodesulfobacteriota_1_draft_adapt_2024_10_12_draft_adapt_validatedbyGapseqVerify_241204_flux_range_rescaled1000to1.RDS',sep=''))

############# User-define the following lines using your models on Katana ############# 

sink(paste0(getwd,'/',new_folder, 'iteration_',auto_num, ".txt"))
print("====================================================================================================")
print(paste("=============================", Sys.time(), '=======================================', sep = ' '))
print(paste("============================= auto_num=", auto_num, '=======================================', sep = ' '))
print("====================================================================================================")
# bacterium_SOB <- BacArena::Bac(my_model,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
arena <- BacArena::Arena(n=grid_no,m=grid_no,Lx=0.012,Ly=0.012) # This grid length represents there are 100 cell distance between two microbial cells.
inocc <- arena@n * arena@m * Inocc_no
bacterium1        <- BacArena::Bac(model1,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium1@type<-"alt_Areni" ######to avoid making bacarena consider our model as the same, this should be finally made by delete part of code in Arena.R
bacterium2        <- BacArena::Bac(model2,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium2@type<-"bin_6"
bacterium3        <- BacArena::Bac(model3,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium3@type<-"bin_16"
bacterium4        <- BacArena::Bac(model4,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium4@type<-"bin_22"
bacterium5        <- BacArena::Bac(model5,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium5@type<-"CCAlphaproteobacteria_3"
bacterium6        <- BacArena::Bac(model6,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium6@type<-"CCAlphaproteobacteria_7"
bacterium7        <- BacArena::Bac(model7,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium7@type<-"CCDehalococcoidia_1"
bacterium8        <- BacArena::Bac(model8,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium8@type<-"CCGammaproteobacteria_3"
bacterium9        <- BacArena::Bac(model9,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium9@type<-"CCGammaproteobacteria_9"
bacterium10        <- BacArena::Bac(model10,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium10@type<-"CCLitorimonas_sp_1"
bacterium11        <- BacArena::Bac(model11,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium11@type<-"CCNitrosopumilus_sp_1"
bacterium12        <- BacArena::Bac(model12,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium12@type<-"CCNitrospirales_2"
bacterium13        <- BacArena::Bac(model13,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium13@type<-"CCPseudohongiella_sp_1"
bacterium14        <- BacArena::Bac(model14,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium14@type<-"CCRhizobiaceae_merged"
bacterium15        <- BacArena::Bac(model15,       lyse = T,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=setAllExInf_value,limit_growth = T)
bacterium15@type<-"CCThermodesulfobacteriota_1"

arena <- BacArena::addOrg(object = arena, specI = bacterium1, amount = inocc*0.0568) # amount can be added by fraction.
arena <- BacArena::addOrg(object = arena, specI = bacterium2, amount = inocc*0.0137)
arena <- BacArena::addOrg(object = arena, specI = bacterium3, amount = inocc*0.0706)
arena <- BacArena::addOrg(object = arena, specI = bacterium4, amount = inocc*0.0118)
arena <- BacArena::addOrg(object = arena, specI = bacterium5, amount = inocc*0.0283)
arena <- BacArena::addOrg(object = arena, specI = bacterium6, amount = inocc*0.0163)
arena <- BacArena::addOrg(object = arena, specI = bacterium7, amount = inocc*0.0353)
arena <- BacArena::addOrg(object = arena, specI = bacterium8, amount = inocc*0.0125)
arena <- BacArena::addOrg(object = arena, specI = bacterium9, amount = inocc*0.0178) # amount can be added by fraction.
arena <- BacArena::addOrg(object = arena, specI = bacterium10, amount = inocc*0.0199)
arena <- BacArena::addOrg(object = arena, specI = bacterium11, amount = inocc*0.2083)
arena <- BacArena::addOrg(object = arena, specI = bacterium12, amount = inocc*0.1046)
arena <- BacArena::addOrg(object = arena, specI = bacterium13, amount = inocc*0.3195)
arena <- BacArena::addOrg(object = arena, specI = bacterium14, amount = inocc*0.0520)
arena <- BacArena::addOrg(object = arena, specI = bacterium15, amount = inocc*0.0326)
# arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange[1], difunc = "pde",
#                            pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM[1]*NutConPercent/100, unit = "mM", add = F,addAnyway = TRUE) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange, difunc = "pde",
                           pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM*NutConPercent/100, unit = "mM", add = F,addAnyway = TRUE) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.

arena@tstep <- time_step
simulation <- BacArena::simEnv(object = arena, time = 1, sec_obj = "opt_ex", continue = T,with_shadow = TRUE)

saveRDS(simulation, file = paste(getwd,'/',new_folder,'/BacArena_',otu_name,'_400grids_mineral_sw_',auto_num,'.RDS',sep = ''))

