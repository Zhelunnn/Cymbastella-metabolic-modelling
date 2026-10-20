################################################################# katana ####################################################################################################
############ 
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

######################## argument ########################
option_list = list(
  optparse::make_option(c("-o", "--script_version"),                 type="character",      default='community_1226_opt_ex',       help="script_version used to distinguish"),
  optparse::make_option(c("-D", "--infile_BacarenaRDS_diet"),  type="character",
                        help="Diet RDS file for bacarena with full path"),
  optparse::make_option(c("-k", "--keywd"),                    type="character",      default="sw_v20210716",    help="(Nutirent/date) keyword used in the folder name of output files, defalt: sw_v20210716"),
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

script_version        = opt$script_version
infile_diet     = opt$infile_BacarenaRDS_diet
keywd           = opt$keywd
grid_no   = opt$arena_mn
death_r   = opt$death_r
Inocc_no  = opt$inocc_no
Cl_no     = opt$cl_no
setAllExInf_value = opt$setAllExInf_value
NutConPercent   = opt$NutConPercent
iter_no   = opt$iter
rm_rate   = opt$rm_rate
time_step = opt$tstep
auto_num  = opt$auto_num
difspeed = opt$difspeed

#####################################################################################################################################################################

diet <- readRDS(infile_diet)
#####################################################################################################################################################################
getwd <- getwd()
getwd

# create a new folder:
new_folder = paste(script_version,'_',keywd,'_death_',death_r,'_rmRate_',rm_rate,'_NutPer_',NutConPercent,'_Inoc_',Inocc_no, '_tstep_',time_step,'h_difspeed_',difspeed,'/',sep = '')
# dir.create(paste(getwd,'/',new_folder,sep = ''))

# auto_num = 2
auto_num_pre = as.character(as.numeric(auto_num) - 1)
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'BacArena_',script_version,'_400grids_mineral_sw_',auto_num_pre,'.RDS',sep = ''))


sink(paste0(getwd,'/',new_folder, 'iteration_',auto_num, ".txt"))
print("====================================================================================================")
print(paste("=============================", Sys.time(), '=======================================', sep = ' '))
print(paste("============================= auto_num=", auto_num, '=======================================', sep = ' '))
print("====================================================================================================")

arena2 <- BacArena::getArena(simulation_loop, 1) # Add the arena in the simulist of 1st iteration.
# for common remove rate (because of water flow)
if (rm_rate!=0){
  arena2@orgdat <- arena2@orgdat[-sample(nrow(arena2@orgdat), round(nrow(arena2@orgdat) * rm_rate/100)), ] #if rm_rate = 10, meaning to remove ~10% of individuals randomly
}


############# ############# User-define the following lines using your models on Katana ############# ############# 
arena2 <- BacArena::addSubs(object = arena2, mediac = diet$Exchange, difunc = "pde",
                            pde = "Diff2d", difspeed = difspeed, smax = diet$Input_mM*NutConPercent/100, unit = "mM", add = F) #Replenish all nutrients by "replacing" (add=F meaning replacing while add=T meaning summing up) after part of population removed.
### difspeed only set the compounds in the medium but keep organisms released compounds to be default

arena2@tstep <- time_step
############# ############# ############# ############# ############# ############# ############# ############# 

simulation <- BacArena::simEnv(object = arena2, time = 1, sec_obj = "opt_ex", continue = T,with_shadow = T) #pFBA for 1 iter

saveRDS(simulation, file = paste(getwd,'/',new_folder,'/BacArena_',script_version,'_400grids_mineral_sw_',auto_num,'.RDS',sep = ''))

