
#' AttributeRiskForRecordI_multiSynthproduct
#'
#' @import matrixStats
#' @export
AttributeRiskForRecordI_multiSynthproduct = function(modelFormulas, i, origdata, syndata, 
                                                     posteriorMCMCs, syntype, H = 50, G = 10, 
                                                     percentBounds = c(0.1, 0.1), 
                                                     additiveBounds = NULL, bounds = NULL,
                                                     guesses = NULL) {
  X_i_syn = list()
  X_i_org = list()
  y_i = list()
  y_i_guesses = list()
  orig_mean = list()
  D = c()
  is_synthesized = c()
  synthesized_predictors = list()
  
  first = TRUE
  
  for (j in 1:length(modelFormulas)) {
    ff = as.formula(modelFormulas[[j]])
    
    utils::str(model <- model.frame(ff, syndata))
    X_i_syn[[j]] <- model.matrix(ff, model)
    utils::str(model <- model.frame(ff, origdata))
    X_i_org[[j]] <- model.matrix(ff, model)
    
    is_synthesized[j] = paste(text = modelFormulas[[j]]$formula[[2]])
    temp = c()
    for (k in 1:length(is_synthesized)) {
      if (is_synthesized[k] %in% paste(text=modelFormulas[[j]]$formula[[3]])) {
        temp = append(temp, is_synthesized[k])
      }
    }
    synthesized_predictors[[j]] = temp
    y_i[[j]] = as.numeric(origdata[i, paste(text = modelFormulas[[j]]$formula[[2]])])
    
    
    
    
    if (j == 2) {
      #return(list(X_i_org[[j]], as.matrix((X_i_org[[j]][i, ]))))
    }
    
    orig_mean[[j]] = t(as.matrix(X_i_org[[j]][i, ])) %*% t(as.matrix(posteriorMCMCs[[j]][, !names(posteriorMCMCs[[j]]) %in% c("sigma")]))
    
    
    
    
    
    
    
    y_i_guesses[[j]] = NULL
    if (paste(text = modelFormulas[[j]]$formula[[2]]) %in% names(Filter(is.factor, origdata))) {
      y_i_guesses[[j]] = as.numeric(levels(origdata[, paste(text = modelFormulas[[j]]$formula[[2]])]))
      D[j] = length(y_i_guesses[[j]])
    } else {
      D[j] = G
      if (is.null(guesses) == FALSE) {
        y_i_guesses[[j]] = guesses
      } else if (is.null(additiveBounds) == FALSE) {
        y_i_guesses[[j]] = seq(y_i[[j]] - additiveBounds[1], y_i[[j]] + additiveBounds[2], length.out = G)
      } else if (is.null(bounds) == FALSE) {
        y_i_guesses[[j]] = seq(bounds[1], bounds[2], length.out = G)
      } else {
        y_i_guesses[[j]] = seq(y_i[[j]] * (1 - percentBounds[1]), y_i[[j]] * (1 + percentBounds[2]), length.out = G)
      }
    }
    
    
    
    if (!(y_i[[j]] %in% y_i_guesses[[j]])) {
      #if (first == TRUE) {
      #  y_i_guesses[[j]] = append(y_i_guesses[[j]], y_i[[j]])
      #  G = G + 1
      #  first = FALSE
      #} else {
      warning("Replaced median value in guess range with true value.")
        y_i_guesses[[j]] = replace(y_i_guesses[[j]], y_i_guesses[[j]] == median(y_i_guesses[[j]]), y_i[[j]])
      #}
    }
  }
  
  CU_i_logZ_all <-rep(NA, prod(D))
  for (j in 1:length(CU_i_logZ_all)) {
    
    currentGuesses = index_to_guesses(j, D, y_i_guesses, is_synthesized)
    #print(j)
    #print(currentGuesses)
    guessed = guessed_mean(X_i_org[[1]][i, ], currentGuesses, synthesized_predictors, posteriorMCMCs, 1, -1)
    
    q_sum_H = (densityCalc(y_i_guesses[[1]][((j-1) %% D[1])+1], syntype[1], list(guessed, posteriorMCMCs[[1]][, "sigma"], D[1]))
               /densityCalc(y_i[[1]], syntype[1], list(orig_mean[[1]], posteriorMCMCs[[1]][, "sigma"], D[1])))
    
    
    #t(as.matrix(X_i_org[[j]][i, ])) %*% t(as.matrix(posteriorMCMCs[[j]][, !names(posteriorMCMCs[[j]]) %in% c("sigma")]))
    
    if (length(modelFormulas) > 1) {  
      for (l in 2:length(modelFormulas)) {
        guessed = guessed_mean(X_i_org[[l]][i, ], currentGuesses, synthesized_predictors, posteriorMCMCs, l, -1)
        
        q_sum_H = q_sum_H * (densityCalc(y_i_guesses[[l]][get_index(D, j, l)], syntype[l], list(guessed, posteriorMCMCs[[l]][, "sigma"], D[l]))
                             /densityCalc(y_i[[l]], syntype[l], list(orig_mean[[l]], posteriorMCMCs[[l]][, "sigma"], D[l])))
        
      }
    }
    
    q_sum_H = sum(q_sum_H)
    log_pq_h_all = rep(NA, H)
    for(h in 1:H) {
      log_p_h = 1
      for (l in 1:length(modelFormulas)) {
        log_p_h = log_p_h * densityCalc(as.numeric(syndata[, paste(text = modelFormulas[[l]]$formula[[2]])]), syntype[l], list(
          as.matrix(X_i_syn[[l]]) %*% t(as.matrix(posteriorMCMCs[[l]][h , !names(posteriorMCMCs[[l]]) %in% c("sigma")])),
          posteriorMCMCs[[l]][h, "sigma"]), D[l])
      }
      log_p_h = sum(log(log_p_h))
      
      guessed = guessed_mean(X_i_org[[1]][i, ], currentGuesses, synthesized_predictors, posteriorMCMCs, 1, h)
      
      log_q_h = (densityCalc(y_i_guesses[[1]][((j-1) %% D[1])+1],syntype[1],list(
                  guessed,
                  posteriorMCMCs[[1]][h, "sigma"], D[1]))
                 /densityCalc(y_i[[1]], syntype[1], list(
                   t(as.matrix(X_i_org[[1]][i, ])) %*% t(as.matrix(posteriorMCMCs[[1]][h, !names(posteriorMCMCs[[1]]) %in% c("sigma")])),
                  posteriorMCMCs[[1]][h, "sigma"], D[1])))
      
      if (length(modelFormulas) > 1) {
        for (l in 2:length(modelFormulas)) {
          guessed = guessed_mean(X_i_org[[l]][i, ], currentGuesses, synthesized_predictors, posteriorMCMCs, l, h)

          log_q_h = log_q_h * (densityCalc(y_i_guesses[[l]][get_index(D, j, l)],syntype[l],list(
                                guessed,
                                posteriorMCMCs[[l]][h, "sigma"], D[l]))
                               /densityCalc(y_i[[l]], syntype[l], list(
                                t(as.matrix(X_i_org[[l]][i, ])) %*% t(as.matrix(posteriorMCMCs[[l]][h, !names(posteriorMCMCs[[l]]) %in% c("sigma")])),
                                posteriorMCMCs[[l]][h, "sigma"], D[l])))
        }
      }
      
      log_q_h = log(log_q_h/q_sum_H)
      
      log_pq_h_all[h] = log_p_h + log_q_h
    }
    
    CU_i_logZ_all[j] = logSumExp(log_pq_h_all)
  }
  
  prob <-exp(CU_i_logZ_all- max(CU_i_logZ_all))/sum(exp(CU_i_logZ_all- max(CU_i_logZ_all)))
  #outcome = array(prob, dim = rep(G, length(modelFormulas)), dimnames = y_i_guesses)
  outcome = array(prob, dim = sapply(y_i_guesses, length), dimnames = y_i_guesses)
  
  marginals = c()
  true_val_string = "outcome["
  if (length(modelFormulas) > 1) {
    for (j in 1:(length(modelFormulas)-1)) {
      true_val_string = paste(true_val_string, "\"", y_i[[j]], "\",", sep = "")
      #marginals[j] = sum(eval(parse(text=paste("outcome[", strrep(",", j - 1), y_i[[j]], strrep(",", length(modelFormulas) - j), "]", sep = ""))))
    }
  }
  true_val_string = paste(true_val_string,"\"", y_i[[length(modelFormulas)]], "\"]", sep = "")
  #marginals[length(modelFormulas)] = sum(eval(parse(text=paste("outcome[", paste(strrep(",", length(modelFormulas) - 1), y_i[[length(modelFormulas)]], "]", sep = "")))))
  
  return(list(FullProb = outcome, TrueMarginals = marginals, TrueValProb = eval(parse(text = true_val_string))))
  
}


densityCalc = function (x, type, otherArgs) {
  if (type == "norm") {
    return(dnorm(x, otherArgs[[1]], otherArgs[[2]]))
  } else if (type == "binom") {
    y = dbinom(x - 1, 1, logistic(otherArgs[[1]]))
    return(y)
  } else if (type == "multinom") {
    return(dmultinom(x, size = otherArgs[[3]], logistic(otherArgs[[1]])))
  } else if (type == "pois") {
    return(dpois(x, exp(otherArgs[[1]])))
  } else {
    stop(paste("Unknown variable type", type))
  }
}

logistic = function (x) {
  return(1/(1+exp(-x)))
}

get_index = function(D, j, l) {
  cnt = 1
  if (l > 1) {
    for (i in 1:(l-1)) {
      cnt = cnt * D[i]
    }
  }
  rtn = as.integer((j - 1) / cnt) + 1
  return(rtn)
}

index_to_guesses = function(j, D, guesses, guess_names) {
  rtn = c()
  for (i in 1:length(D)) {
    if (i == 1) {
      rtn[i] = guesses[[i]][((j-1) %% D[i]) + 1] 
    } else {
      rtn[i] = guesses[[i]][as.integer((j-1) / D[i]) + 1]
    }
  }
  names(rtn) = guess_names
  return(rtn)
}

guessed_mean = function(X_row, currentGuesses, synthesized_predictors, posteriorMCMCs, l, h) {
  temp = X_row
  for (k in 1:length(synthesized_predictors[[1]])) {
    temp = replace(temp, synthesized_predictors[[1]][k], currentGuesses[synthesized_predictors[[1]][k]])
  }
  if (h == -1) {
    guessed = t(as.matrix(temp)) %*% t(as.matrix(posteriorMCMCs[[l]][, !names(posteriorMCMCs[[l]]) %in% c("sigma")]))
  } else {
    guessed = t(as.matrix(temp)) %*% t(as.matrix(posteriorMCMCs[[l]][h, !names(posteriorMCMCs[[l]]) %in% c("sigma")]))
  }
  return(guessed)
}
