### polyDFE v1.0: predicting DFE and alpha from polymorphism data
### Copyright (c) 2016  Paula Tataru
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
### 
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
### Contact: paula@birc.au.dk 


source("postprocessing.R")

########################### 
# parse the polyDFE output
########################### 
est = parseOutput("output/example_5_basin")
print(est)

########################### 
## estimate alpha
########################### 
# from DFE
estimateAlpha(est[[1]])
# from the DFE, but use the S_adv limit like in Galtier 2016 paper
estimateAlpha(est[[1]], supLimit=2)
# use divergence counts
# for that, I need to read the divergence data
div = parseDivergenceData('input/example_5')
print(div)
estimateAlpha(est[[1]], div=div)
# now use the S_adv limit
estimateAlpha(est[[1]], div=div, supLimit=2)
estimateAlpha(est[[1]], div=div, supLimit=5)
# account for the misattributed polymorphism
estimateAlpha(est[[1]], div=div, poly=T)
# and with S_adv
estimateAlpha(est[[1]], div=div, poly=T, supLimit=5)

########################### 
# LRT and AIC
########################### 
# this works both by comparing directly parsed estimates
# or by comparing two files
est1 = parseOutput("output/example_5_")
est2 = parseOutput("output/example_5_joint")
# it compares run i in est1 with run i in est2
compareModels(est1, est2)
compareModels('output/example_1_basin', 'output/example_1_aut')

########################### 
# making init lines
########################### 
# sometimes we want the best found parameters from an optimization run
# to be used as the starting parameters for a new run
# for that, we need to create the right init file needed for polyDFE
# eps_cont is a parameter that models contamination of neutral sites with selected sites
# it's not identifiable - so it, in fact, should never be estimated
createInitLines('output/example_1_aut', 'input/test', 
					 startingID = 1, fix = c("eps_cont", "a"))
createInitLines('output/example_5_basin', 'input/test', 
					 startingID = 2, fix = c("eps_cont", "a"))
# now automatically group the r parameters
createInitLines('output/example_1_aut', 'input/test', 
					 startingID = 3, fix = c("eps_cont", "a"), 
					 groupingDiff = 0.001)
createInitLines('output/example_5_basin', 'input/test', 
					 startingID = 4, fix = c("eps_cont", "a"), 
					 groupingDiff = 0.001)
# the test_grouping file is empty - no r parameters where close enough
# try a larger distance
createInitLines('output/example_1_aut', 'input/test', 
					 startingID = 5, fix = c("eps_cont", "a"), 
					 groupingDiff = 0.1)
createInitLines('output/example_5_basin', 'input/test', 
					 startingID = 6, fix = c("eps_cont", "a", "S_p -10"), 
					 groupingDiff = 0.1)
# the ids in test_init and test_grouping are made to match
# when grouping r values that are already grouped, the groups obtained in
# test_grouping are not correct, so should not be used
# the R code doesn't know how the r where originally grouped 
