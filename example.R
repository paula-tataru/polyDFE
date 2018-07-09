### polyDFE v2.0: predicting DFE and alpha from polymorphism data
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

# multiple runs of polyDFE can be stored in the same list
# for easier processing
est = list()
polyDFEfiles = c("output/example_2_del", "output/example_2_full_A",
                 "output/example_2_nonuis_C", "output/example_2_nonuis_s_C",
                 "output/example_2_nonuis_sj_C", "output/example_2_novar_C",
                 "output/example_2_init_C", "output/example_2_full_C")
for (filename in polyDFEfiles)
{
    est = c(est, parseOutput(filename))
}
print(length(est))

# what are the gradients?
grad = sapply(est, function(e) e$criteria)
print(grad)
# for which files are the gradients a bit large?
polyDFEfiles[which(grad > 0.01)]
# for those runs it might help to run some basin hopping iterations on top

# on what input files was polyDFE ran?
print(sapply(est, function(e) e$input))

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
colnames(dfe) = toNames(c(-100, -10, -1, 0, 1))
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
print(estimateAlpha(est[[8]]))
# alpha_dfe, but use the S_adv limit like in Galtier 2016
print(estimateAlpha(est[[8]], supLimit = 2))
# alpha_div
# for that, divergence data has to be read first
div = parseDivergenceData('input/example_2')
print(div)
print(estimateAlpha(est[[8]], div = div))
# alpha_div with S_adv limit
print(estimateAlpha(est[[8]], div = div, supLimit = 2))
# a higher supLimit reduces the value of alpha
print(estimateAlpha(est[[8]], div = div, supLimit = 5))
# turn off correction for misattributed polymorphism
print(estimateAlpha(est[[8]], div = div, poly = FALSE))
print(estimateAlpha(est[[8]], div = div, poly = FALSE, supLimit = 5))
# calculate alpha for all runs of polyDFE found in est
# and compare with the values returned from polyDFE
# there will be small differences 
# because one is calculated in C and the other in R
alpha = sapply(est, function(e) c("polyDFE alpha_dfe" = unname(unlist(e$alpha)["alpha_dfe"]),
                                  "R alpha_dfe" = estimateAlpha(e),
                                  "polyDFE alpha_div" = unname(unlist(e$alpha)["alpha_div"]),
                                  "R alpha_div" = estimateAlpha(e, div = div)))
colnames(alpha) = basename(polyDFEfiles)
print(alpha)

###########################
# model testing
###########################
# this works both by comparing two files
print(compareModels("output/example_2_del", "output/example_2_full_C"))
# or directly parsed estimates
print(compareModels(est[1], est[8]))
# if multiple runs are given, it compares run i in est1 with run i in est2
# for example, we can compare the del DFE with the full DFEs
print(compareModels(est[c(1, 1)], est[c(2, 8)]))
# by default, the function detects that runs 1 and 2 are not nested
# as run 1 was a deleterious DFE obtained from model C
# but this is, in fact, nested in model A too
# we can then enforce nestedness
print(compareModels(est[c(1, 1)], est[c(2, 8)], nested = TRUE))
# if only est1 is given, it just returns the AIC
aic = compareModels(est)
print(aic)
# "output/example_2_nonuis_C" gives the best AIC
# which also fits best with the simulations, where the r vlues where 1
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
# createInitLines is not aware of the ID's already present in the file
# so startingID should be set carefully to make sure
# the IDs present in the file are unique
# eps_cont is a parameter that models contamination of neutral sites with selected sites
# it's not identifiable - so createInitLines always write its flag to be 1
createInitLines('output/example_2_full_C', 'input/test', startingID = 1)
# the function can be given a list of polyDFE runs
# but they have to have the same DFE model
createInitLines(est, 'input/test', startingID = 2)
createInitLines(est[3:8], 'input/test', startingID = 2)
# init files can have parameters for different models
createInitLines('output/example_2_full_A', 'input/test', startingID = 8)
# fix (flag 1) mutation variability
createInitLines(est[3:8], 'input/test', startingID = 9, fix = c("eps_cont", "a"))
# fix all parameters - this is useful for just calculating the likelihood
# for a specific test of parameters, sometimes needed in model testing
createInitLines(est[3:8], 'input/test', startingID = 15, fix = "all")
# we can automatically group the r parameters
createInitLines(est[8], 'input/test', startingID = 21,
                groupingDiff = 0.001)
# no neighbouring r values are closer than 0.001, so the test_grouping file is not updated
# we can try a larger distance
createInitLines(est[8], 'input/test', startingID = 22,
                groupingDiff = 0.1)


#################################################################################
# model testing for sharing parameters
#################################################################################
est = c(parseOutput("output/example_1_2_indep"),
        parseOutput("output/example_1_2_share"))
# what are the gradients?
grad = sapply(est, function(e) e$criteria)
print(grad)
# they are a bit big, running basin hopping might help

# on what input files was polyDFE ran?
print(sapply(est, function(e) e$input))
# each entry in est has 2 inputs as polyDFE was ran on 2 files jointly

# what is the discretized DFE?
getDiscretizedDFE(est[[1]])
getDiscretizedDFE(est[[2]])
# just using sapply as before doesn't work anymore
# as it doesn't return the right format
# because getDiscretizedDFE now returns a matrix
# containg 2 rows, one for each input file that was used
# we can use lapply and rbind instead the result
do.call("rbind", lapply(est, getDiscretizedDFE))
# the last 2 rows are the same because those result form the 
# polyDFe run where the DFE was shared by the 2 datasets

# what is alpha?
est[[1]]$alpha
# I have, again, 2 sets of alphas, one for each input dataset
do.call("rbind", est[[1]]$alpha)
do.call("rbind", 
        lapply(est, function(e) do.call("rbind", e$alpha)))
# the last 2 rows of alpha_dfe are the same because those result form the 
# polyDFe run where the DFE was shared by the 2 datasets
# when calling estimateAlpha, it returns 2 estimates, 
# one for each input dataset
estimateAlpha(est[[1]])
# divergence then has to be supplied as a list to estimate alpha_div
div = lapply(est[[1]]$input, parseDivergenceData)
estimateAlpha(est[[1]], div = div)
# calculate alpha for all runs of polyDFE found in est
# and compare with the values returned from polyDFE
alpha = lapply(est, 
               function(e) 
                   rbind("polyDFE alpha_dfe" = do.call("rbind", e$alpha)[, "alpha_dfe"],
                         "R alpha_dfe" = estimateAlpha(e),
                         "polyDFE alpha_div" = do.call("rbind", e$alpha)[, "alpha_div"],
                         "R alpha_div" = estimateAlpha(e, div = div)))
alpha = do.call("rbind", alpha)
print(alpha)

# test if the DFE should be shared or not
compareModels(est[1], est[2])
# p-value says that we should use different DFEs!

# when creating init lines, two line are written
# one for each intput dataset
createInitLines(est[1], 'input/test', startingID = 23,
                groupingDiff = 0.1)
