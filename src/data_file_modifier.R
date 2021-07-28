# STILL NOT FINISHED
# data_file_modifier.R
# Date created:  October 31, 2015
# This code reads in ctl and dat files for modifying

# Store the file names from which data is available
filename <- vector( length=5)
filename[1]="PWS_ASA.dat"
filename[2]="PWS_ASA(ESS).ctl"
filename[3]="PWS_ASA(age_comps).ctl"
filename[4]="PWS_ASA(surveys).ctl"
filename[5]="PWS_ASA(disease).ctl"

# Here is my function to read in the data from each file, and the headers (separately) as a list
setwd("~/PWS_herring/Model/ASA_recent_Data/TMB_2015")
source(file="~/PWS_herring/Model/ASA_recent_Data/TMB_2015/data_reader.R")
source(file="~/PWS_herring/Model/ASA_recent_Data/TMB_2015/data_header_reader.R")

#
setwd("~/PWS_herring/Model/ASA_recent_Data/extended_Age0")
PWS_ASA.dat <- data_reader(filename=filename[1])
headers <- data_header_reader(filename=filename[1])

# START HERE
new.wd <- paste0(modelPath,"output_plots")
dir.create(path=new.wd)
newdir <- "~/PWS_herring/Model/ASA_recent_Data/extended_Age0"

nyr <- PWS_ASA.dat[[1]]
indexing <- 1
for(j in 1:length){
  data <- PWS_ASA.dat[[j]]
  stop.here <- which(headers[[2]])[j]
  bits <- headers[[1]][indexing:stop.here]
  if(j==1){
    write(bits,file=paste0(newdir,filename[1])),append=FALSE)
    write(data,file=paste0(newdir,filename[1])),sep="\t",append=TRUE)
  }else {
    write(x,file=paste0(newdir,filename[1])),sep="\t",append=TRUE)
  }
  
  
}


Data = list("w_a_a"=data_reader(filename=filename[1])[[3]], # From PWS_ASA.dat
            "fecun"=data_reader(filename=filename[1])[[4]],
            "pc"=data_reader(filename=filename[1])[[5]],
            "pk"=data_reader(filename=filename[1])[[6]],
            "fbc"=data_reader(filename=filename[1])[[7]],
            "gc"=data_reader(filename=filename[1])[[8]],
            "sc"=data_reader(filename=filename[1])[[9]],
            "f_sp"=data_reader(filename=filename[1])[[10]],
            "Z_3_8"=data_reader(filename=filename[1])[[11]],
            # From PWS_ASA(ESS).ctl
            "ESS_Se"=data_reader(filename=filename[2])[[1]], 
            "ESS_Sp"=data_reader(filename=filename[2])[[2]],
            # From PWS_ASA(age_comps).ctl
            "seine"=data_reader(filename=filename[3])[[1]], 
            "spac"=data_reader(filename=filename[3])[[2]],
            # From PWS_ASA(surveys).ctl
            "mdm"=data_reader(filename=filename[4])[[1]], 
            "egg"=data_reader(filename=filename[4])[[2]],
            "cv_egg"=data_reader(filename=filename[4])[[3]],
            "hydADFG_start"=data_reader(filename=filename[4])[[4]],
            "hydADFG"=data_reader(filename=filename[4])[[5]],
            "hydPWSSC_start"=data_reader(filename=filename[4])[[6]],
            "hydPWSSC"=data_reader(filename=filename[4])[[7]],
            "cv_hydPWSSC"=data_reader(filename=filename[4])[[8]],
            # From PWS_ASA(disease).ctl
            "vhsv"=data_reader(filename=filename[5])[[5]], 
            "ich"=data_reader(filename=filename[5])[[6]])

