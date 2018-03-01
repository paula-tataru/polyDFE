### polyDFE v1.0: predicting DFE and alpha from polymorphism data
### Copyright (c) 2018  Paula Tataru
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
# parse the output from one run of polyDFE
est = parseOutput("output/example_2_full_C")
# est is a list of length one
print(length(est))
# where entry 1 contains the information from the polyDFE run
print(names(est[[1]]))
# which parameters where estimated?
print(est[[1]]$estimated)
print(est[[1]])

# setting init = TRUE forces the function to return
# the values that were estimated using starting values
# given by the user with -i
est_1 = parseOutput("output/example_2_init_C")
est_2 = parseOutput("output/example_2_init_C", init = TRUE)
cat(est_1[[1]]$lk, "\t", est_1[[1]]$grad, "\t", est_1[[1]]$values["S_d"], "\n")
cat(est_2[[1]]$lk, "\t", est_2[[1]]$grad, "\t", est_2[[1]]$values["S_d"], "\n")

# multiple runs of polyDFE can be stored in the same list
# for easier processing
est = list()
polyDFEfiles = c("output/example_2_del", "output/example_2_full_A",
                 "output/example_2_nonuis_C", "output/example_2_novar_C",
                 "output/example_2_init_C", "output/example_2_full_C")
for (filename in polyDFEfiles)
{
    est = c(est, parseOutput(filename))
}
print(length(est))
# what are the gradients?
print(sapply(est, function(e) e$grad))
# none of them are particullary large, so the optimization was most likely succesfull

########################### 
## summarize the DFE
########################### 
getDiscretizedDFE(est[[2]])
# change the binning
getDiscretizedDFE(est[[2]], sRanges = c(-100, -10, -1, 0))
# discretize the DFE for all runs of polyDFE found in est
dfe = t(sapply(est, getDiscretizedDFE))
# we can also change the binning
dfe = t(sapply(est, getDiscretizedDFE, sRanges = c(-100, -10, -1, 0, 1)))
# we can give names to the dfe matrix that partialy reflect which model was ran
rownames(dfe) = sapply(est, getModelName)
print(dfe)
# the discretized DFE can also be visualized as a barplot
barplot(dfe, beside = TRUE, legend.text = TRUE)
# using a log scale for the y-axis esnables an easier visual comparison
# for that, we have to replace the 0 values in dfe with NA
dfe[which(dfe == 0)] = NA
barplot(dfe, beside = TRUE, log = 'y', col = 1:nrow(dfe))
# we can use the original filenames for the legend instead
legend("topright", legend = basename(polyDFEfiles), fill = 1:nrow(dfe))

########################### 
## estimate alpha
########################### 
# alpha_dfe
print(estimateAlpha(est[[6]]))
# alpha_dfe, but use the S_adv limit like in Galtier 2016
print(estimateAlpha(est[[6]], supLimit = 2))
# alpha_div
# for that, divergence data has to be read first
div = parseDivergenceData('input/example_2')
print(div)
print(estimateAlpha(est[[6]], div = div))
# alpha_div with S_adv limit
print(estimateAlpha(est[[6]], div = div, supLimit = 2))
# a higher supLimit reduces the value of alpha
print(estimateAlpha(est[[6]], div = div, supLimit = 5))
# turn off correction for misattributed polymorphism
print(estimateAlpha(est[[6]], div = div, poly = FALSE))
print(estimateAlpha(est[[6]], div = div, poly = FALSE, supLimit = 5))
# calculate alpha for all runs of polyDFE found in est
# alpha_dfe
print(sapply(est, estimateAlpha))
# alpha_div
print(sapply(est, estimateAlpha, div = div))

########################### 
# model testing
########################### 
# this works both by comparing tow files
print(compareModels("output/example_2_del", "output/example_2_full_C"))
# or directly parsed estimates
print(compareModels(est[1], est[6]))
# it multiple runs are given, it compares run i in est1 with run i in est2
# for example, we can compare the del DFE with the full DFEs
print(compareModels(est[c(1, 1)], est[c(2, 6)]))
# by default, the function detects that runs 1 and 2 are not nested
# as run 1 was a deleterious DFE obtained from model C
# but this is, in fact, nested in model A too
# we can then enfornce nestedness
print(compareModels(est[c(1, 1)], est[c(2, 6)], nested = TRUE))
# if only est1 is given, it just returns the AIC
aic = compareModels(est)
print(aic)
# "output/example_2_nonuis_C" gives the best AIC
# which also fits with the simulations, where the r vlues where 1
print(polyDFEfiles[which.min(aic$AIC[, "AIC"])])

########################### 
# model averaging
########################### 
# calculate the AIC weights
aic_weights = getAICweights(est)
print(aic_weights)
# calculate model-averaged DFE
avg_dfe = rowSums(sapply(1:length(est), 
                         function(i) aic_weights[i, "weight"] * dfe[i, ]), na.rm = TRUE)
print(avg_dfe)
barplot(avg_dfe)

########################### 
# making init lines
########################### 
# sometimes we want the best found parameters from an optimization run
# to be used as the starting parameters for a new run
# for that, we need to create the right init file for polyDFE
# eps_cont is a parameter that models contamination of neutral sites with selected sites
# it's not identifiable - so it should never be estimated
createInitLines('output/example_2_full_C', 'input/test', startingID = 1)
# the function can be given a list of polyDFE runs
# but they have to have the same DFE model
createInitLines(est, 'input/test', startingID = 2)
createInitLines(est[3:5], 'input/test', startingID = 2)
# init files can have parameters for different models
createInitLines('output/example_2_full_A', 'input/test', startingID = 5)
# fix (flag 1) mutation variability
createInitLines(est[3:6], 'input/test', startingID = 6, fix = c("eps_cont", "a"))
# fix all parameters - this is useful for just calculating the likelihood 
# for a specific test of parameters, sometimes needed in model testing
createInitLines(est[3:6], 'input/test', startingID = 10, fix = "all")
# we can automatically group the r parameters
createInitLines(est[6], 'input/test', startingID = 14,
                groupingDiff = 0.001)
# no neighbouring r values are closer than 0.001, so the test_grouping file is not updated
# we can try a larger distance
createInitLines(est[6], 'input/test', startingID = 15,
                groupingDiff = 0.1)
# when grouping r values that are already grouped, the groups obtained in
# test_grouping are not correct, so should not be used
# the R code doesn't know the original grouping
