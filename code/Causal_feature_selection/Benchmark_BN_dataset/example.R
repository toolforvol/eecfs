
#  http://www.bnlearn.com/bnrepository/


# ************ Names of the 31 networks are as follows ************ 

# ----- Discrete Bayesian Networks

# Small Networks (<20 nodes)   asia, cancer, earthquake, sachs, survey

# Medium Networks (20?C50 nodes)   alarm, barley, child, insurance, mildew, water

# Large Networks (50?C100 nodes)   hailfinder, hepar2, win95pts

# Very Large Networks (100?C1000 nodes)   andes, diabetes, link, pathfinder, pigs, munin1, munin2, munin3, munin4';

# Massive Networks (>1000 nodes)   munin

# ----- Continues Bayesian Networks (Gaussian)

# Medium Networks (20?C50 nodes) ecoli70, magic-niab

# Large Networks (50?C100 nodes) magic-irri

# Very Large Networks (101?C1000 nodes) arth150

# ----- Continues Bayesian Networks (Conditional Linear Gaussian)

# Small Networks (<20 nodes) sangiovese

# Medium Networks (20?C50 nodes) mehra-original, mehra-complete



cat("\014")  
rm(list = ls())

source('Gen_Data.R')

data_name<-"link"

data_samples<-5000      

Gen_Data(data_name,data_samples)

# data and graph will be generated in the data folder



