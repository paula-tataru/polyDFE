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

mergeNames = function(x, y)
{
	pos = which(y == x, arr.ind = TRUE)
	if (length(pos) == 0)
	{
		return(y)
	}
	newY = NULL
	if (min(pos) > 1)
	{
		# copy everything from before the first occurence of x
		newY = y[1:(min(pos) - 1)]
	}
	# merge the names by merging the entries where x is found
	# with the neighbour entries to the right
	newY = c(newY, paste(y[pos], y[pos + 1], sep = ' '))
	if (max(pos) + 2 < length(y))
	{
		# copy everthing from after the last occurence of x
		newY = c(newY, y[(max(pos) + 2):length(y)])
	}
	return(newY)
}

splitBySpace = function(x)
{
	# replace all \t with space
	y = gsub("\t", " ", x)
	# split by space
	y = strsplit(y, split = ' ')[[1]]
	# remove empty entries
	y = y[which(y != '', arr.ind = TRUE)]
	# remove '--' if it's the first entry
	if (y[1] == '--')
	{
		y = y[-1]
	}
	# if I have r or S_p parameters in here, they are separated by spaces
	# merge them back
	y = mergeNames('r', y)
	y = mergeNames('S_p', y)
	return(y)
}

findPos = function(x, lines)
{
	which(unlist(lapply(lines, grepl, pattern = x, fixed = TRUE)),
	      arr.ind = TRUE)
}

strip = function(s)
{
	# strip leading and trailing white spaces
	gsub("^\\s+|\\s+$", "", s)
}

removeEmpty = function(x)
{
	# remove from x the entries that are empty
	x[which(unlist(lapply(strip(x), nchar)) > 0, arr.ind = TRUE)]
}

splitEqual = function(x)
{
	strsplit(x, split = " = ", fixed = T)[[1]]
}

extractInfo = function(lines, init = FALSE)
{
    # remove empty lines
    lines = removeEmpty(lines)
    
    values = NULL
    
    # find out on what file I ran polyDFE
    pos = findPos('-- Performing inference', lines)[1]
    filename = strsplit(lines[pos], split = ' ')[[1]][5]
    
    # find out which are the estimated parameters
    # the info is in the line following -- Starting local optimization
    pos = findPos('-- Starting local optimization', lines)[1] + 1
    
    if (is.na(pos))
    {
        # I didn't have any -- Starting local optimization
        # nothing was estimated then!
        # look for ---- Joint likelihood
        pos = findPos('---- Joint likelihood', lines)
        # extract the likelihood
        lk = strsplit(lines[pos], split = ' ')[[1]]
        lk = as.numeric(lk[length(lk)])
        estimated = NA
        grad = NA
        # extract the model, it's in position 2
        model = strsplit(lines[pos + 1], split = ': ')[[1]][2]
    } else
    {
        # figure out what parameters have been estimated
        estimated = splitBySpace(lines[pos])
        # the first and last 4 contain
        # iteration number, ln lk, gradient and status
        estimated = estimated[2:(length(estimated) - 4)]
        
        # find out what is the best likelihood found
        # what model I used
        # and what are the parameters
        # get the position of best likelihood
        pos = findPos('---- Best joint likelihood', lines)
        
        # extract the likelihood and gradient
        lk = strsplit(lines[pos], split = ' ')[[1]]
        # remove empty entries
        lk = lk[which(unlist(lapply(lk, nchar)) > 0, arr.ind = TRUE)]
        grad = as.numeric(lk[length(lk)])
        lk = as.numeric(lk[length(lk) - 3])
        # extract the model, it's in position 2
        model = strsplit(lines[pos + 1], split = ': ')[[1]][2]
    }
    
    # extract the rest of the parameters; I have 6 lines left
    for (i in seq(from = pos + 2, to = pos + 6, by = 2))
    {
        # get the parameters on this line
        params = c(names(values), splitBySpace(lines[i]))
        # get the values from the next line
        values = c(values, as.numeric(splitBySpace(lines[i + 1])))
        names(values) = params
    }
    
    if (init)
    {
        # I actually want the parameters estimated from the init values
        pos = findPos('---- Using provided initial values', lines)
        if (length(pos) == 0)
        {
            # no init was used
            warning("polyDFE was not run with the provided initial values")
            return(NULL)
        }
        # move on to starting the optimization
        # using the initial values is always the first one
        pos = findPos('-- Starting local optimization', lines)[1] + 1
        if (!is.na(pos))
        {
            # find the last line of the output
            # it has to end with a line that starts with --
            end = findPos('--', lines)
            end = end[which(end > pos)][1]
            # find the last line that had only numbers
            for (i in end:pos)
            {
                aux = strsplit(lines[i], split = ' ')[[1]]
                # get rid of empty entries
                aux = removeEmpty(aux)
                aux = tryCatch(as.numeric(aux),
                               warning = function(w) return(w))
                if (typeof(aux) != "list")
                {
                    # valid conversion, get results
                    # the first and last 3 contain
                    # iteration number, lnlk, gradient and status
                    values[estimated] = aux[2:(length(aux) - 3)]
                    # update likelihood and gradient
                    lk = as.numeric(aux[length(aux) - 2])
                    grad = as.numeric(aux[length(aux) - 1])
                    break
                }
            }
        } 
    }
    
    # add eps_an and eps_cont if missing
    if (!"eps_an" %in% names(values))
    {
        values = c("eps_an" = 0, values)
    }
    if (!"eps_cont" %in% names(values))
    {
        values = c(values[1], "eps_cont" = 0, values[-1])
    }
    # add p_b and S_b if missing
    # only for models B - C, which also contain paramter b
    if (!"p_b" %in% names(values) && (model == 'B' || model == 'C'))
    {
        pos = which(names(values) == "b")
        if (length(pos) > 0)
        {
            values = c(values[1:pos],
                       "p_b" = 0,
                       "S_b" = 0,
                       values[(pos + 1):length(values)])
        }
    }
    # add a if missing
    if (!"a" %in% names(values))
    {
        pos = which(names(values) == "theta_bar")
        values = c(values[1:pos],
                   "a" = -1,
                   values[(pos + 1):length(values)])
    }
    # add lambda if missing
    if (!"lambda" %in% names(values))
    {
        pos = which(names(values) == "eps_cont")
        values = c(values[1:pos],
                   "lambda" = 0,
                   values[(pos + 1):length(values)])
    }
    
    # extract expectations
    # expectations are in the order
    # P_neut, P_sel, D_neut, D_sel, mis_neut, mis_sel
    pos = findPos('---- Expected', lines)
    # first, figure out n
    n = pos[2] - pos[1]
    # also, figure out where expectations end
    pos = c(pos[1], pos[2], pos[2] + n)
    expec = NULL
    for (i in 1:2)
    {
        expec = rbind(expec,
                      unlist(lapply(lines[(pos[i] + 1):(pos[i + 1] - 1)],
                                    function(x)
                                        as.numeric(splitEqual(x)[2]))))
    }
    expec = cbind(expec, c(NA, NA), c(NA, NA))
    colnames(expec) = c(paste0("E[P_z(", 1:(n - 1), ")]"), "E[D_z]", "E[mis_z]")
    rownames(expec) = c("neut", "sel")
    # add E[D_neut] and E[D_sel]
    pos = findPos('E[D_neut] = ', lines)
    if (length(pos) == 1)
    {
        expec[1, n] = as.numeric(splitEqual(lines[pos])[2])
        expec[2, n] = as.numeric(splitEqual(lines[pos + 1])[2])
    }
    # add E[mis_neut] and E[mis_sel]
    pos = findPos('E[mis_neut] = ', lines)
    if (length(pos) == 1)
    {
        expec[1, n + 1] = as.numeric(splitEqual(lines[pos])[2])
        expec[2, n + 1] = as.numeric(splitEqual(lines[pos + 1])[2])
    }
    
    # extract alpha
    pos = findPos('---- alpha', lines)
    alpha = NULL
    for (p in pos)
    {
        aux = splitEqual(lines[p])
        alpha = c(alpha, as.numeric(aux[2]))
        names(alpha)[length(alpha)] = strsplit(aux[1], split = "---- ")[[1]][2]
    }
    
    return(
        list(
            input = filename,
            model = model,
            lk = lk,
            grad = grad,
            values = values,
            estimated = estimated,
            n = n,
            expec = expec,
            alpha = alpha
        )
    )
}

parseOutput = function(filename, init = FALSE)
{
	# extract the model, estimated parameters and likelihood
	# from the file
	# it can potentially contain multiple runs
	# so find all of them
	# a run starts wtih -- Starting local optimization
	# and ends with --- Best joint likelihood
	# I need to look at the starting of the optimization
	# to see which paramters are actually estimated
	
	f = file(filename)
	lines = readLines(f)
	close(f)
	
	# find which lines contain ---- Performing inference on
	pos = findPos('---- Performing inference on', lines)
	
	# add last line at the end
	pos = c(pos, length(lines))
	
	# go over each run and extract the info
	res = lapply(1:(length(pos) - 1),
					 function(i) extractInfo(lines[pos[i]:(pos[i + 1] - 1)], init))
	
	return(res)
}

parseSFSData = function(filename)
{
	f = file(filename)
	lines = readLines(f)
	close(f)
	
	# figure out how many lines to skip - i.e. where I have the i j n info
	# these are comment lines starting with # or empty lines
	# when checking for empty lines,
	# I need to strip leading and trailing white spaces
	i = 1
	while (nchar(strip(lines[i])) == 0 || substr(lines[i], 1, 1) == "#")
	{
		i = i + 1
	}
	
	data = read.table(filename, skip = i, comment.char = "#")
	
	info_line = na.omit(as.numeric(unlist(strsplit(lines[i], "[^0-9]+"))))
	
	# sum the fragments
	data = rbind(colSums(data[1:info_line[1],]),
					 colSums(data[(info_line[1] + 1):(info_line[1] + info_line[2]),]))
	
	# give it names
	rownames(data) = c("neut", "sel")
	sfs = paste0("p(", 1:(info_line[3] - 1), ")")
	if (info_line[3] + 2 == ncol(data))
	{
		# with divergence data
		colnames(data) = c(sfs, "length_sfs", "d", "length_d")
	} else
	{
		# without divergence data
		colnames(data) = c(sfs, "length_sfs")
	}
	
	return(data)
}

parseDivergenceData = function(filename)
{
	data = parseSFSData(filename)
	
	if (!("length_d" %in% colnames(data)))
	{
		warning("File does not contain divergence data")
		return(NA)
	}
	
	# now get divergence counts
	div = c("length_neut" = data["neut", "length_d"], 
			  "length_sel" = data["sel", "length_d"], 
			  "obs_neut" = data["neut", "d"], 
			  "obs_sel" = data["sel", "d"])
	
	return(div)
}

getQuantiles = function(est)
{
    model = est$model
    values = est$values
    
    if (model == "D")
    {
        return(list(q = NULL, p = NULL, m = NULL))
    }
    
    sep = 0.0005
    p = seq(from = 0.0, to = 1, by = sep)
    p[1] = 1 - 0.99999999
    p[length(p)] = 0.99999999
    n = length(p)
    q = c()
    # I have to store the corresponding probabilities too
    prob = c()
    if (model != 'D')
    {
        if (model == 'A')
        {
            S_max = values["S_max"]
            S_bar = values["S_bar"]
        } else
        {
            S_max = 0
            S_bar = values["S_d"]
        }
        q = S_max - qgamma(p, shape = values["b"], 
                           scale = (S_max - S_bar) / values["b"], 
                           lower.tail = FALSE)
        # last quantile should in fact be S_max
        q[length(q)] = S_max
        prob = rep(sep, n - 1)
        if (model == 'B')
        {
            q = c(q, values['S_b'])
            prob = prob * (1 - values["p_b"])
            prob = c(prob, values["p_b"])
        } else if (model == 'C')
        {
            q = c(q[2:length(q)], qexp(p, rate = 1 / values['S_b']))
            prob = prob * (1 - values["p_b"])
            prob = c(prob, rep(sep, n - 1) * values["p_b"])
        }
    }
    # make the quantiles unique
    newQ = q[1]
    newP = prob[1]
    for (i in 2:(length(q) - 1))
    {
        if (q[i] == q[i - 1])
        {
            newP[length(newP)] = newP[length(newP)] + prob[i]
        } else
        {
            newQ = c(newQ, q[i])
            newP = c(newP, prob[i])
        }
    }
    # deal with the last interval - increase the last quantile if equal
    i = length(q)
    if (q[i] == q[i - 1])
    {
        q[i] = q[i] + 1e-7
    } 
    newQ = c(newQ, q[i])
    q = newQ
    prob = newP
    
    # check probabilities sum to 1
    if (!isTRUE(all.equal(sum(prob), 1)))
    {
        stop("Numerical integration failed, probabilities don't sum up to 1.")
    }
    
    # calculate mid-points
    mid = (q[1:(length(q) - 1)] + q[2:length(q)]) / 2
    
    return(list(q = q, p = prob, m = mid))
}

altIntegrate = function(f, lower, upper, quantiles, ...)
{
    q = quantiles$q
    # sum up the probabilties for the quantiles that are between lower and upper
    pos = which(lower <= q & q <= upper)
    # mid and prob i is for interval [i, i + 1]
    # that means that I always have to skip the last entry in pos
    pos = pos[-length(pos)]
    value = sum(quantiles$p[pos] * f(quantiles$m[pos], ...))
    # now I have to deal with the cases where lower and upper
    # are found between quantiles and I need to split an interval
    if (lower > q[1] && !(lower %in% q))
    {
        i = which(lower > q)
        i = length(i)
        prop = (q[i + 1] - lower) / (q[i + 1] - q[i])
        value = value + prop * quantiles$p[i] * f(quantiles$m[i], ...)
    }
    if (upper < q[length(q)] && !(upper %in% q))
    {
        i = which(upper < q)
        i = i[1] - 1
        prop = (upper - q[i]) / (q[i + 1] - q[i])
        value = value + prop * quantiles$p[i] * f(quantiles$m[i], ...)
    }
    
    return(value)
}

extractSP = function(values)
{
	sp = unlist(strsplit(names(values[grep('S_p', names(values))]), split = ' '))
	sp = as.numeric(sp[seq(from = 2, to = length(sp), by = 2)])
	return(sp)
}

toNames = function(ranges)
{
	names = paste('<', ranges[1])
	for (i in 1:(length(ranges) - 1))
	{
		names = c(names, paste('(', ranges[i], ', ', ranges[i + 1], ')', sep=''))
	}
	names = c(names, paste(ranges[length(ranges)], '<'))
	return(names)
}

cdfRefDisGamma = function(S, S_bar, b, S_max = 0)
{
    pgamma(S_max - S, shape = b, scale = (S_max - S_bar) / b, 
           lower.tail = FALSE)
}

getDiscretizedDFE = function(estimates, sRanges = c(-100, -10, -1, 0, 1, 10))
{
    # if the list doesn't contain model
    # but it's entry does
    # throw an error
    if (!("model" %in% names(estimates)) && ("model" %in% names(estimates[[1]])))
    {
        stop("estimates should contain just one entry ",
             "from the list returned by parseOutput")
        return(NULL)
    }
    # add -Inf and Inf to range so I can calculate integrals
    # make sure sRanges is sorted
    sRanges = sort(sRanges)
    sRanges = c(-Inf, sRanges, Inf)
    
    model = estimates$model
    values = estimates$values
    
    # the get the discretized DFE, I use the CDF functions
    if (model == "A")
    {
        # displaced reflected Gamma
        cdf = cdfRefDisGamma(sRanges, values["S_bar"], values["b"], values["S_max"])
    } else if (model == "B" || model == "C")
    {
        # I have to add 0 to sRanges if it's not part of it
        newRanges = sort(unique(c(sRanges, 0)))
        pos = which(newRanges == 0)
        # displaced Gamma for the negative side
        cdfDel = ((1 - values["p_b"]) 
                  * cdfRefDisGamma(newRanges[1:pos], values["S_d"], values["b"]))
        # positive side
        if (model == "B")
        {
            cdfBen = rep(values["p_b"], length(newRanges) - pos)
            pos = which(newRanges[pos:length(newRanges)] <= values["S_b"])
            pos = pos[length(pos)] - 1
            cdfBen[1:pos] = 0
        } else
        {
            # exponenetial for the positive side
            cdfBen = (values["p_b"] 
                      * pexp(newRanges[(pos + 1):length(newRanges)], 
                             rate = 1 / values["S_b"]))
        }
        cdf = c(cdfDel, cdfDel[length(cdfDel)] + cdfBen)
        # if 0 was not part of the ranges, I have to fix that
        if (!(0 %in% sRanges))
        {
            cdf = c(# the deleterious part up to the value before 0
                cdfDel[1:(length(cdfDel) - 1)], 
                # the beneficial part 
                cdfDel[length(cdfDel)] + cdfBen)
        }
    } else
    {
        # model D
        cdf = rep(0, length(sRanges))
        # extract the S_p values
        sp = extractSP(values)
        for (i in 1:length(sp))
        {
            pos = which(sRanges > sp[i])
            cdf[pos] = cdf[pos] + values[paste("S_p", sp[i])]
        }
    }
    dfe = cdf[2:length(cdf)] - cdf[1:(length(cdf) - 1)]
    
    # check it sums to almost 1
    if (!isTRUE(all.equal(sum(dfe), 1)))
    {
        warning("Discretized DFE sums to ", sum(dfe), 
                "\n\tDiscretized DFE has been rescaled to sum to 1.")
    }
    dfe = dfe / sum(dfe)
    
    names(dfe) = toNames(sRanges[2:(length(sRanges) - 1)])
    return(dfe)
}

getModelName = function(estimates)
{
	modelName = estimates$model
	if (estimates$model == "A")
	{
		if (!"S_max" %in% estimates$estimated && estimates$values["S_max"] == 0)
		{
			modelName = paste(modelName, "del")
		}
	}
	if (estimates$model == "B" || estimates$model == "C")
	{
		if (!"p_b" %in% estimates$estimated && estimates$values["p_b"] == 0)
		{
			modelName = paste(modelName, "del")
		}
	}
	pos = grep("r ", names(estimates$values))
	if (!"r" %in% estimates$estimated 
		 	&& isTRUE(all.equal(estimates$values[pos], 
		 	                    rep(1, length(pos)), 
		 	                    check.attributes = FALSE)))
	{
		modelName = paste(modelName, "no r", sep = ", ")
	}
	if (!"eps_an" %in% estimates$estimated && estimates$values["eps_an"] == 0)
	{
		modelName = paste(modelName, "no eps", sep = ", ")
	}
	
	return(modelName)
}

getDivergenceExpectation = function(S)
{
    # this is conditional on a given S
    pos = which(S != 0)
    expec = rep(1, length(S))
    if (length(pos) > 0)
    {
        expec[pos] = S[pos] / (1 - exp(-S[pos]))
    }
    # when using my own numerical integration, I do not multiply with the density
    return(expec)
}

getMisattributedPolyExpectation = function(S, n)
{
    # this is conditional on a given S
    # I replace the integral over x with the Taylor expansion
    # just like in polyDFE
    pos0 = which(S == 0)
    if (length(pos0) != 0)
    {
        newS = S[-pos0]
    } else
    {
        newS = S
    }
    pos = which(newS >= 0)
    if (length(newS) > 0)
    {
        # try numerical integration
        aux = sapply(newS,
                     function(s)
                         tryCatch(
                             integrate(function(x)
                                 x^(n - 1) * (1 - exp(-s * (1 - x))) / (1 - x),
                                 0, 1)$value,
                             warning = function(war)
                             {
                                 return(NA)
                             },
                             error = function(err)
                             {
                                 return(NA)
                             }))
        # for the positions where numerical integration failed
        # I use the Taylor approximation
        failed = which(is.na(aux))
        if (length(failed) > 0)
        {
            # Taylor approximation
            term = rep(1, length(newS[failed]))
            aux[failed] = rep(1, length(newS[failed]))
            for (k in 0:10)
            {
                term = term * (-newS[failed]) / (n + k)
                aux[failed] = aux[failed] + term / (k + 1)
            }
            aux[failed] = - aux[failed]
        }
        if (length(pos) == 0)
        {
            aux = aux * exp(newS) / (exp(newS) - 1)
        } else if (length(pos) == length(newS))
        {
            aux = aux / (1 - exp(-newS))
        } else
        {
            aux[pos] = aux[pos] / (1 - exp(-newS[pos]))
            aux[-pos] = aux[-pos] * exp(newS[-pos]) / (exp(newS[-pos]) - 1)    
        }
    } else
    {
        aux = NULL
    }
    if (length(pos0) != 0)
    {
        # for S == 0, I just use the neutral
        # re-instert it at the original position of 0
        aux = append(aux, 1 / n, after = pos0 - 1)
    }
    # when using my own numerical integration, I do not multiply with the density
    return(aux)
}

hasPosSel = function(est)
{
    model = est$model
    values = est$values
	if (model == 'A')
	{
		return(!(is.na(values['S_max'] || values['S_max'] == 0)))
	}
	if (model == 'B' || model == 'C')
	{
		return(!(is.na(values['p_b']) || is.na(values['S_b'])
					|| values['p_b'] == 0 || values['S_b'] == 0))
	}
	if (model == 'D')
	{
		return(any(extractSP(values) > 0))
	}
}

estimateAlpha = function(estimates, supLimit = 0, div = NULL, poly = TRUE)
{
    # if the list doesn't contain model
    # but it's entry does
    # throw an error
    if (!("model" %in% names(estimates)) && ("model" %in% names(estimates[[1]])))
    {
        stop("estimates should contain just one entry ",
             "from the list returned by parseOutput")
        return(NULL)
    }
    
    model = estimates$model
    values = estimates$values
    n = estimates$n
    
    # make sure supLimit is >= 0
    if (supLimit < 0)
    {
        warning("supLimit cannot be negative, setting it to 0")
        supLimit = 0
    }
    
    # I need to use the critical values for a better integration
    quantiles = getQuantiles(estimates)
    
    # make sure I have p_b and S_b, they might be missing when not estimated
    if (model != 'D')
    {
        if (!'p_b' %in% names(values))
        {
            values = c(values, 0, 0)
            names(values) = c(names(values)[1:(length(values) - 2)], 
                              "p_b", "S_b")
        }	
    }
    
    # calculate the expected alpha from the estimated DFE
    # it's dependent on the model
    # for the integral giving the divergence expectations
    if (model != 'D')
    {
        divDel = altIntegrate(getDivergenceExpectation,
                              -Inf, supLimit, quantiles)
    } else
    {
        # extract the S_p values
        sp = extractSP(values)
        posSupLim = which(sp <= supLimit)
        if (length(posSupLim) > 1)
        {
            posSupLim = posSupLim[length(posSupLim)]
        } else
        {
            posSupLim = length(sp)
        }
        # calculate the deletrious divergence
        divDel = 0
        for (i in 1:posSupLim)
        {
            divDel = (divDel 
                      + getDivergenceExpectation(sp[i]) 
                      * values[paste("S_p", sp[i])])
        }
    }
    
    if (is.null(div))
    {
        # the whole integral should be multiplied by lambda
        # but it cancels out in the division at the end
        # if no positive selection, return o
        if (!hasPosSel(estimates))
        {
            return(0)
        }
        
        # for models B and D, I do not have an integral
        if (model != 'B' && model != 'D')
        {
            divBen = altIntegrate(getDivergenceExpectation,
                                  supLimit, Inf, quantiles) 
        } else
        {
            if (model == 'B')
            {
                divBen = getDivergenceExpectation(values['S_b']) * values["p_b"]	
            } else
            {
                divBen = 0
                if (posSupLim < length(sp))
                {
                    for (i in (posSupLim + 1):length(sp))
                    {
                        divBen = (divBen + getDivergenceExpectation(sp[i]) 
                                  * values[paste("S_p", sp[i])])
                    }	
                } else
                {
                    return(0)
                }
            }
        }
        
        return(divBen / (divDel + divBen))
    }
    
    # if I have divergence data,
    # I calcualte alpha according to observed divergence counts
    # if poly = T, I account for the misattributed polymorphism
    # for this, I need to have theta and r_n
    if (poly)
    {
        # if lambda was not estimated, I do not have r_n, so just set it to 1
        if (!"lambda" %in% names(values))
        {
            r_n = 1
        } else
        {
            # figure out what is my biggest r
            pos = grep("r ", names(values), fixed = T)
            r_n = values[pos[length(pos)]]
        }
        
        # get the misattributed neutral polymorphism
        fake_div_neut = div["length_neut"] * values["theta_bar"] * r_n / n
        # get the misattributed selected polymorphism
        if (model != 'D')
        {
            fake_div_sel = altIntegrate(getMisattributedPolyExpectation,
                                        -Inf, Inf, quantiles, n)
        } else
        {
            fake_div_sel = 0
            for (i in 1:length(sp))
            {
                fake_div_sel = (fake_div_sel 
                                + getMisattributedPolyExpectation(sp[i], n) 
                                * values[paste("S_p", sp[i])])
            }
        }
        
        fake_div_sel = div["length_sel"] * values["theta_bar"] * r_n * fake_div_sel
        
        div["obs_neut"] = div["obs_neut"] - fake_div_neut
        div["obs_sel"] = div["obs_sel"] - fake_div_sel
    }
    alpha = as.numeric(1 - 
                       ((div["length_sel"] / div["length_neut"])
                        * (div["obs_neut"] / div["obs_sel"])
                        * divDel))
    
    return(alpha)
}

getAIC = function(estimates)
{
	return(2 * length(estimates$estimated) - 2 * estimates$lk)
}

compareModels = function(est1, est2 = NULL, nested = NULL)
{
	# compare the two sets of estimated parameters
	# get the AIC and, if the models are nested,
	# calculate the LRT
	# if the files contain multiple entries
	# I compare across coresponding entries
	
	if (is.character(est1))
	{
		est1 = parseOutput(est1)
	}
	if (!is.null(est2) & is.character(est2))
	{
		est2 = parseOutput(est2)
	}
	
	if (!is.null(est2) & length(est1) != length(est2))
	{
		warning('Estimates do not contain the same number of runs: ',
			 length(est1), ' vs ', length(est2), '\n')
	}
	
	mm = length(est1)
	if (!is.null(est2))
	{
		mm = max(mm, length(est2))
	}
	
	aic = NULL
	lrt = NULL
	for (i in 1:mm)
	{
		aic1 = NA
		aic2 = NA
		if (i <= length(est1))
		{
			aic1 = getAIC(est1[[i]])
		}
		
		if (is.null(est2))
		{
			aic = rbind(aic, c(length(est1[[i]]$estimated), 
									 est1[[i]]$lk, aic1))
			next
		}
		
		if (i <= length(est2))
		{
			aic2 = getAIC(est2[[i]])
		}
		
		aic = rbind(aic, c(length(est1[[i]]$estimated), 
								 est1[[i]]$lk, aic1, 
								 length(est2[[i]]$estimated), 
								 est2[[i]]$lk, aic2))
		
		if ((i > length(est1)) || i > length(est2))
		{
			break
		}
		
		# check if models are nested
		# if the r's do not match, then this will return FALSE
		# but then I can control that from the "outside" with the nested option
		# models should be the same
		nest = (est1[[i]]$model == est2[[i]]$model)
		# estimated parameters are nested
		nest = (nest 
				  && (all(est1[[i]]$estimated %in% est2[[i]]$estimated)
				  	   || all(est2[[i]]$estimated %in% est1[[i]]$estimated)))
		# fixed parameters should be equal
		fixed = intersect(setdiff(names(est1[[i]]$values), est1[[i]]$estimated),
								setdiff(names(est2[[i]]$values), est2[[i]]$estimated))
		nest = (nest && isTRUE(all.equal(est1[[i]]$values[fixed],
													est2[[i]]$values[fixed])))
		# if models have the same number of parameters, 
		# they can't be nested
		nest = (nest && 
				  	length(est1[[i]]$estimated) != length(est2[[i]]$estimated))
		
		# calculate LRT
		# skip estimates where the number of parameters is the same
		# these cannot possibly be nested
		if (est1[[i]]$input != est2[[i]]$input)
		{
		    warning("Input files are different for run", i)
		    lrt = rbind(lrt, c(NA, est1[[i]]$lk, est2[[i]]$lk, NA))
		} else if (isTRUE(nested) || (nest && is.null(nested)))
		{
		    # the D statistic is positive
		    # if I don't care which is the bigger model,
		    # I just take the absolute value
		    D = 2 * abs(est1[[i]]$lk - est2[[i]]$lk)
		    # D should follow a chi-square distribution
		    # with degrees of freedom given by the number of extra parameters
		    # estimated in the larger model
		    df = abs(length(est1[[i]]$estimated) - length(est2[[i]]$estimated))
		    lrt = rbind(lrt, 
		                c(df, est1[[i]]$lk, est2[[i]]$lk, 
		                  pchisq(D, df, lower.tail = FALSE)))
		} else
		{
			lrt = rbind(lrt, c(NA, est1[[i]]$lk, est2[[i]]$lk, NA))
		}
	}
	
	if (!is.null(est2))
	{
		colnames(aic) = c('df model 1', 'log lk model 1', 'AIC model 1', 
								'df model 2', 'log lk model 2', 'AIC model 2')

		if (is.null(nested) || isTRUE(nested))
		{
			colnames(lrt) = c('df', 'log lk model 1', 
									'log lk model 2', 'p-value')	
			return(list(AIC = aic, LRT = lrt))
		}
		return(list(AIC = aic))
	}
	
	colnames(aic) = c('df', 'log lk', 'AIC')
	return(list(AIC = aic))
}

getAICweights = function(estimates)
{
	aic = compareModels(estimates)$AIC
	# calculate Delta AIC
	aic[, "AIC"] = aic[, "AIC"] - min(aic[, "AIC"])
	# now calculate the weights
	w = exp(-0.5 * aic[, "AIC"])
	aic = cbind(aic, "weight" = w / sum(w))
	colnames(aic)[3] = "delta AIC"
	return(aic)
}

grouping = function(rValues, diff = NA)
{
    if (is.na(diff))
    {
        # no grouping
        return(list(groups = names(rValues), r = rValues))
    }
    
	# calculate grouping of r values
	# by merging values that have an absolute difference of at most diff
	groups = NULL
	i = 1
	j = 2
	while (i <= length(rValues))
	{
		while (j <= length(rValues) && abs(rValues[i] - rValues[j]) <= diff)
		{
			j = j + 1
		}
		# from i to j-1, I have a group
		# however, the indexes are off by 1, so I write j instead of j-1
		groups = c(groups, j)
		# update i and j
		i = j
		j = j + 1
	}
	# now calculate the new r values by taking the mean from each group
	newValues = NULL
	# add the start of the groups - which is r2
	auxGroups = c(1, groups)
	for (i in 1:length(groups))
	{
		# re-adjust the indexes
		pos = auxGroups[i]:(auxGroups[i + 1] - 1)
		newValues = c(newValues, round(mean(rValues[pos]), digits = 5))
	}
	names(newValues) = paste("r", groups)
	return(list(groups = groups, r = newValues))
}

writeModelD = function(values, rPos, fixedInfo = NULL)
{
	# create proper init lines for model D
	# it is different than the rest of the models
	pos = grep('S_p', names(values))
	pos = pos[c(1, length(pos))]
	if (is.null(fixedInfo))
	{
		# use the names
		s = sprintf('%-6s', paste0('S', 1:(pos[2] - pos[1] + 1)))
		p = paste0('p', 1:(pos[2] - pos[1] + 1))
		auxSP = paste(sprintf('%-16s', paste(s, p)), collapse = ' ')
		toWrite = paste(sprintf('%-12s', names(values)[1:(pos[1] - 1)]),
							 collapse = ' ')
		toWrite = paste(toWrite, auxSP, 'r', '\n')
		return(toWrite)
	}
	
	# extract the S_p values
	sp = extractSP(values)
	
	toWrite = paste(fixedInfo[1:(pos[1] - 1)],
						 sprintf('%-10g', values[1:(pos[1] - 1)]),
						 collapse = ' ')
	toWrite = paste(toWrite, 
						 paste(fixedInfo[pos[1]:pos[2]],
						 		sprintf('%-4g', sp),
						 		sprintf('%-9g', values[pos[1]:pos[2]]),
						 		collapse = ' '))
	
	return(toWrite)
}

createInitLines = function(estimates, outputfile, startingID = 1, 
                           fix = c("eps_cont"),
                           groupingDiff = NA)
{
	# create init lines with the estimated parameters from inputfile or estimates
	# using ids starting at the given value
	# and fix the given parameters
	# write to outputfile
	
	if (is.character(estimates))
	{
		allEst = parseOutput(estimates)
	} else
	{
		allEst = estimates
	}
	
	id = startingID
	model = allEst[[1]]$model
	# check all runs in allEst use the same model
	for (est in allEst)
	{
	    if (est$model != model)
	    {
	        stop('Not all polyDFE runs use the same model!\n')
	    }
	}
	allParams = names(allEst[[1]]$values)
	
	# update fix if it says "all"
	if (fix[1] == "all")
	{
	    fix = allParams
	}
	
	# I need to treat the r parameters a bit differently
	# as I have just one indicator for them!
	rPos = which("r " == substr(allParams, 1, 2))[1]
	# the parameters are already in the correct order
	fixedInfo = rep(0, rPos)
	# update fix first - if it contains r
	pos = which(fix == 'r')
	fix[pos] = allParams[grep('r ', allParams)[1]]
	aux = which(unlist(lapply(allParams[1:rPos], function(x) x %in% fix)))
	fixedInfo[aux] = rep(1, length(aux))
	
	# before writing all the parameters
	# write an ID line
	initFile = paste(outputfile, "init", sep = "_")
	cat(sprintf('%-6s', '# ID'), file = initFile, append = TRUE)
	
	if (model != 'D')
	{
		toWrite = paste(sprintf('%-12s', allParams[1:(rPos - 1)]),
							 collapse = ' ')
		toWrite = paste(toWrite, 'r', '\n')	
	} else
	{
		toWrite = writeModelD(allEst[[1]]$values, rPos)
	}
	cat(toWrite, file = initFile, append = TRUE)
	
	writeGroupHeader = FALSE
	
	for (est in allEst)
	{
		if (model != 'D')
		{
			toWrite = paste(fixedInfo[1:(rPos - 1)],
								 sprintf('%-10g', est$values[1:(rPos - 1)]),
								 collapse = ' ')	
		} else
		{
			toWrite = writeModelD(est$values, rPos, fixedInfo)
		}
		
		# calculate the grouping
		myRValues = est$values[rPos:(length(est$values))]
		myGroups = grouping(myRValues, groupingDiff)
		# add the r parametesr
		toWrite = paste(toWrite, paste(c(fixedInfo[rPos], myGroups$r),
												 collapse = ' '))
		# add the id
		toWrite = paste(sprintf('%-5d', id), toWrite, '\n', sep = ' ')
		cat(toWrite, file = initFile, append = TRUE)
		# if I groupped anything or I have to write all groups,
		# write to the grouping file
		if (length(myGroups$groups) < length(myRValues))
		{
			if (!writeGroupHeader)
			{
				groupingFile = paste(outputfile, "grouping", sep = "_")
				cat(paste("# ID no-groups   i j ... k",
							 "# groups are: [r_1] (automated group, as r_1 = 1), ",
							 "# [r_2, r_i], [r_i+1, r_j], ..., [r_k+1, r_{n-1}]",
							 "# where n is the number of sampled sequences\n",
							 sep = "\n"),
					 file = groupingFile, append = TRUE)
				writeGroupHeader = TRUE
			}
			
			# first, ID
			toWrite = sprintf('%-5d', id)
			# then the number of groups
			toWrite = paste(toWrite, sprintf('%-5d', length(myGroups$groups)))
			# and lastly the groups
			toWrite = paste(toWrite, paste(myGroups$groups, collapse = ' '), '\n\n')
			cat(toWrite, file = groupingFile, append = TRUE)
		}
		id = id + 1
	}
	cat('\n', file = initFile, append = TRUE)
}

bootstrapData = function(inputfile, outputfile = NULL, rep = 1)
{
	# this only works for data that has just two fragments
	# one neutral and one selected
	data = parseSFSData(inputfile)
	n = which(colnames(data) == "length_sfs")
	# extract and remove lengths
	l = data[, n, drop = FALSE]
	data = data[, -n]
	if (ncol(data) >= n)
	{
		l = cbind(l, data[, ncol(data), drop = FALSE])
		data = data[, -ncol(data)]
	}
	
	# bootstrap the data
	neut = sapply(data["neut", ], function(lam) rpois(max(rep, 2), lambda = lam))
	sel = sapply(data["sel", ], function(lam) rpois(max(rep, 2), lambda = lam))
	
	# write it to file
	if (is.null(outputfile))
	{
		outputfile = inputfile
	}
	dig = floor(log10(rep))
	for (i in 1:rep)
	{
		filename = paste0(outputfile, "_", formatC(i, digits = dig, flag = "0"))
		cat("# Bootstraped data from", inputfile, "\n", 
			 file = filename, append = FALSE)
		cat("1 1", n, "\n", file = filename, append = TRUE)
		
		# the neutral sfs
		for (z in neut[i, 1:(n - 1)])
		{
			cat(sprintf("%6s", z), " ", file = filename, append = TRUE)
		}
		# the neutral length
		cat(sprintf("%6s", l["neut", "length_sfs"]), 
			 file = filename, append = TRUE)
		# the neutral divergence count and length
		if (ncol(l) == 2)
		{
			cat(sprintf("%6s", neut[i, "d"]), sprintf("%6s", l["neut", "length_d"]), 
				 file = filename, append = TRUE)
		}
		cat("\n", file = filename, append = TRUE)
		
		# the selected sfs
		for (z in sel[i, 1:(n - 1)])
		{
			cat(sprintf("%6s", z), " ", file = filename, append = TRUE)
		}
		# the selected length
		cat(sprintf("%6s", l["sel", "length_sfs"]), 
			 file = filename, append = TRUE)
		# the selected divergence count and length
		if (ncol(l) == 2)
		{
			cat(sprintf("%6s", sel[i, "d"]), sprintf("%6s", l["sel", "length_d"]), 
				 file = filename, append = TRUE)
		}
		cat("\n", file = filename, append = TRUE)
	}
}
