# =============================================================================
# functions.R
#
# Helper functions sourced by the analysis scripts. Packages are loaded by the
# calling scripts: lme4 (Multilmer), nlme and doBy (multigls), FactoMineR
# (GroupingSummary), mixOmics and dplyr (sPLS functions) and stringr
# (compare_metab).
#
#   Pre-processing   Multilmer()
#   MWAS             multigls(), RunLinear()
#   Clustering       GroupingSummary()
#   sPLS             tune_spls_cor(), perf_cor(), moe(), pls.perf.par(),
#                    quality.var()
#   Validation       compare_metab()
# =============================================================================


# Pre-processing ----------------------------------------------------------

#' Linear mixed models fitted to each column of a data matrix
#'
#' Fits `y ~ 1 + fixed + (1 | random)` for each column of `ydata`. Random
#' effects with a variance below `tol` are dropped. Used to remove plate
#' effects ("denoising"): the residuals returned are the observations minus
#' the random intercepts.
#'
#' @param ydata matrix of outcomes (samples x features).
#' @param xdata data frame with the fixed and random effect variables.
#' @param fixed names of the fixed effect variables (optional).
#' @param random names of the random intercept variables.
#' @param var name of a fixed effect to test (continuous or binary only).
#' @param tol tolerance below which a random effect variance is set to zero.
#' @param progressbar show a progress bar.
#' @param residuals return the denoised data.
#' @return list with the number of observations, the coefficients for `var`
#'   (if provided), the random intercepts, the denoised data (if
#'   `residuals = TRUE`) and the intra-class correlation of each random effect.
Multilmer <- function(ydata, xdata, fixed = NULL, random, var = NULL, tol = 1e-5, progressbar = TRUE, residuals = FALSE) {
  if (!is.null(var)) {
    coefficients <- matrix(NA, nrow = ncol(ydata), ncol = 3)
    rownames(coefficients) <- colnames(ydata)
    colnames(coefficients) <- c("coef", "coef.se", "pval")
  }
  nobs <- rep(NA, ncol(ydata))
  names(nobs) <- colnames(ydata)

  if (residuals) {
    Residuals <- matrix(NA, nrow = nrow(ydata), ncol = ncol(ydata))
    rownames(Residuals) <- rownames(ydata)
    colnames(Residuals) <- colnames(ydata)
  }

  # Random intercepts: one matrix (levels x features) per random effect
  random_intercepts <- list()
  for (r in 1:length(random)) {
    mylevels <- unique(as.character(xdata[, random[r]]))
    ri <- matrix(NA, ncol = ncol(ydata), nrow = length(mylevels))
    rownames(ri) <- mylevels
    colnames(ri) <- colnames(ydata)
    random_intercepts[[r]] <- ri
  }

  if (progressbar) {
    progbar <- txtProgressBar(style = 3)
  }

  # Intra-class correlation
  icc <- matrix(NA, nrow = ncol(ydata), ncol = length(random))
  rownames(icc) <- colnames(ydata)
  colnames(icc) <- random

  for (k in 1:ncol(ydata)) {
    if (progressbar) {
      setTxtProgressBar(progbar, k / ncol(ydata))
    }
    x <- ydata[, k]
    nobs[k] <- nrow(na.exclude(cbind(x, xdata[c(fixed, random)])))
    if (!is.null(fixed)) {
      fe <- paste(fixed, collapse = "+")
    }
    re <- paste(paste0("(1|", random, ")"), collapse = "+")
    if (!is.null(fixed)) {
      f <- paste0("x~1+", fe, "+", re)
    } else {
      f <- paste0("x~1+", re)
    }
    model <- lmer(as.formula(f),
      data = xdata, REML = FALSE,
      control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
    )
    vcor <- data.frame(VarCorr(model))
    rownames(vcor) <- vcor[, 1]
    vcov.sum <- sum(vcor[, "vcov"])
    vcor <- vcor[random, "vcov"]
    icc[k, ] <- vcor / vcov.sum

    newrandom <- random[vcor > tol]
    if (length(newrandom) == 0) {
      # No random effect left: linear model
      if (!is.null(fixed)) {
        fe <- paste(fixed, collapse = "+")
        f <- paste0("x~1+", fe)
        model <- lm(as.formula(f), data = xdata)
        if (residuals) {
          resids <- residuals(model) + data.matrix(xdata[, fixed]) %*% coef(model)[fixed] + coef(model)[1]
          Residuals[, k] <- resids
        }
      } else {
        if (residuals) {
          Residuals[, k] <- x
        }
      }
      if (!is.null(var)) {
        newfixed <- fixed[!fixed %in% var]
        if (length(newfixed) == 0) {
          fe <- 1
        } else {
          fe <- paste(newfixed, collapse = "+")
        }
        f <- paste0("x~1+", fe)
        model0 <- lm(as.formula(f), data = xdata)
        myanova <- anova(model, model0)
        coefficients[k, 3] <- myanova$`Pr(>F)`[2]
      }
    } else {
      if (any(vcor < tol)) {
        # Some but not all random effect terms are smaller than tolerance
        re <- paste(paste0("(1|", newrandom, ")"), collapse = "+")
        if (!is.null(fixed)) {
          f <- paste0("x~1+", fe, "+", re)
        } else {
          f <- paste0("x~1+", re)
        }
        model <- lmer(as.formula(f),
          data = xdata, REML = FALSE,
          control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
        )
      }
      for (r in 1:length(newrandom)) {
        absolute_rid <- which(random == newrandom[r])
        rintercepts <- ranef(model)[[newrandom[r]]]
        random_intercepts[[absolute_rid]][rownames(rintercepts), k] <- unlist(rintercepts)
      }
      if (!is.null(var)) {
        newfixed <- fixed[!fixed %in% var]
        fe <- paste(newfixed, collapse = "+")
        re <- paste(paste0("(1|", newrandom, ")"), collapse = "+")
        f <- paste0("x~1+", fe, "+", re)
        model0 <- lmer(as.formula(f),
          data = xdata, REML = FALSE,
          control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
        )
        myanova <- anova(model, model0)
        if (sum(grepl(var, rownames(summary(model)$coefficients))) > 1) {
          stop("This code only works for continuous or binary 'var'. Please create dummy variables.")
        }
        coefficients[k, 3] <- myanova$`Pr(>Chisq)`[2]
      }
      if (residuals) {
        resids <- residuals(model) + data.matrix(xdata[, fixed]) %*% fixef(model)[fixed] + fixef(model)[1]
        Residuals[, k] <- resids
      }
    }
    if (!is.null(var)) {
      coefficients[k, 1:2] <- summary(model)$coefficients[grep(var, rownames(summary(model)$coefficients)), 1:2]
    }
  }
  names(random_intercepts) <- random

  out <- list(nobs = nobs)
  if (!is.null(var)) {
    out$coefficients <- as.data.frame(coefficients, stringsAsFactors = FALSE)
  }
  out$random_intercepts <- random_intercepts
  if (residuals) {
    out$residuals <- Residuals
  }
  out$icc <- icc
  return(out)
}


# MWAS --------------------------------------------------------------------

#' Multivariate normal (MVN) models for repeated measures
#'
#' For each feature, fits a generalised least squares model with an
#' unstructured within-subject correlation and wave-specific variances:
#' `feature ~ Expo.mean + Expo.delta + covars`, where `Expo.mean` is the
#' subject's mean exposure and `Expo.delta` the deviation from it (exposure is
#' first centred on its mean before swimming). Returns the estimate for
#' `Expo.delta`.
#'
#' @param ydata matrix of outcomes (samples x features).
#' @param xdata data frame with columns `ID`, `wave`, `Expo` and `covars`,
#'   in the same row order as `ydata`.
#' @param covars names of the adjustment covariates.
#' @param progressbar not used.
#' @param standardise_outcome divide each feature by its standard deviation.
#' @param complete_case fit the models on complete cases only.
#' @return data frame with one row per feature: `OMIC_ID`, `coef`, `coef.se`,
#'   `pval`, `tval`.
multigls <- function(ydata, xdata, covars = NULL, progressbar = FALSE, standardise_outcome = FALSE, complete_case = FALSE) {
  coefficients <- matrix(NA, nrow = ncol(ydata), ncol = 5)
  colnames(coefficients) <- c("OMIC_ID", "coef", "coef.se", "pval", "tval")

  if (progressbar) {
    progbar <- txtProgressBar(style = 3)
  }

  data <- xdata
  data$Expo <- data$Expo - mean(subset(data, wave == 1)$Expo, na.rm = TRUE) # centre exposure at the mean level before swimming

  # Mean exposure per subject and within-subject deviation
  tmp <- summaryBy(Expo ~ ID, data = data)
  data <- merge(data, tmp, by = "ID", all.x = TRUE)
  data$Expo.delta <- data$Expo - data$Expo.mean

  for (CurrOmicIndex in 1:ncol(ydata)) {
    data$Omic <- ydata[, CurrOmicIndex]
    if (standardise_outcome) {
      data$Omic <- data$Omic / sd(data$Omic)
    }
    gls_formula <- as.formula(paste0("I(Omic) ~ (Expo.mean + Expo.delta + ", paste(covars, collapse = "+"), ")"))

    if (complete_case) {
      model <- gls(gls_formula,
        data = data[complete.cases(data), ], # removing samples with any missing information
        correlation = corSymm(form = ~ 1 | ID),
        weights = varIdent(form = ~ 1 | wave),
        na.action = na.pass
      )
    } else {
      model <- gls(gls_formula,
        data = data,
        correlation = corSymm(form = ~ 1 | ID),
        weights = varIdent(form = ~ 1 | wave),
        na.action = na.pass
      )
    }

    coefficients[CurrOmicIndex, 1] <- colnames(ydata)[CurrOmicIndex]
    coefficients[CurrOmicIndex, 2] <- summary(model)$tTable["Expo.delta", "Value"]
    coefficients[CurrOmicIndex, 3] <- summary(model)$tTable["Expo.delta", "Std.Error"]
    coefficients[CurrOmicIndex, 4] <- summary(model)$tTable["Expo.delta", "p-value"]
    coefficients[CurrOmicIndex, 5] <- summary(model)$tTable["Expo.delta", "t-value"]
  }
  coefficients <- as.data.frame(coefficients)
  coefficients[, 2:5] <- apply(coefficients[, 2:5], 2, as.numeric)
  return(coefficients)
}

#' Linear regression of each outcome on a predictor
#'
#' @param covars data frame with the predictor and confounders.
#' @param predictor name of the predictor (exposure).
#' @param outcomes matrix of outcomes (samples x features).
#' @param conf names of the confounders.
#' @return data frame with one row per outcome: estimate, SE, p-value and 95%
#'   confidence interval for the predictor, and the number of observations.
RunLinear <- function(covars, predictor, outcomes, conf) {
  mysummary.list <- list()
  for (k in 1:ncol(outcomes)) {
    mydatauniv <- cbind(covars[, c(predictor, conf)], y = outcomes[, k])
    mydatauniv <- na.exclude(mydatauniv)
    f <- paste0("y~", paste(c(predictor, conf), collapse = "+"))
    mymodel <- lm(as.formula(f), data = mydatauniv)
    mycoef <- summary(mymodel)$coefficients
    mycoef <- as.data.frame(matrix(mycoef[c(predictor), c(1, 2, 4)], ncol = 3))
    colnames(mycoef) <- c("estimate", "SE", "p-value")
    mysummary <- cbind(
      outcome = colnames(outcomes)[k],
      nobs = nrow(mydatauniv),
      mycoef,
      lci = mycoef[, "estimate"] - 1.96 * mycoef[, "SE"],
      uci = mycoef[, "estimate"] + 1.96 * mycoef[, "SE"]
    )
    rownames(mysummary) <- NULL
    mysummary.list[[k]] <- mysummary
  }
  summary <- do.call(rbind, mysummary.list)
  return(summary)
}


# Clustering --------------------------------------------------------------

#' Summarise groups of features by their first principal component
#'
#' @param x data matrix (samples x features).
#' @param group named vector of group memberships (names are features).
#' @param summary "centroid" (PC1 scores) or "medoid" (feature contributing
#'   most to PC1).
#' @param verbose show a progress bar.
#' @return list with `xsummarised` (samples x groups, scaled) and
#'   `description` (one row per feature: group, pairwise correlations within
#'   the group, variance explained by PC1, loadings and contributions).
GroupingSummary <- function(x, group, summary = "centroid", verbose = TRUE) {
  # Getting data transpose
  x <- t(x)

  # Number of groups
  k <- max(group)

  # Initialisation of the objects
  xsummarised <- matrix(NA, nrow = k, ncol = ncol(x))
  myrownames <- rep(NA, length = k)
  mytable <- data.frame(
    feature = names(group), group = unname(group),
    min_cor = rep(1, length(group)),
    max_cor = rep(1, length(group)),
    mean_cor = rep(1, length(group)),
    ev = rep(NA, length(group)),
    loadings = rep(NA, length(group)),
    contribution = rep(NA, length(group))
  )
  rownames(mytable) <- mytable[, 1]

  # Loop over the groups
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  for (i in 1:k) {
    # Extracting data in the group
    xcluster <- x[names(group)[which(group == i)], , drop = FALSE]

    # PCA on the group
    mypca <- FactoMineR::PCA(t(xcluster), graph = FALSE, ncp = 1)

    # Proportion of explained variance by PC1
    mytable[mytable$group == i, "ev"] <- mypca$eig[1, 2] / 100

    # Loadings coefficients on PC1 (sign chosen so that most loadings are positive)
    if ((sum(mypca$var$coord > 0) / length(mypca$var$coord)) > 0.5) {
      mytable[rownames(xcluster), "loadings"] <- mypca$var$coord / sqrt(mypca$eig[1, 1])
    } else {
      mytable[rownames(xcluster), "loadings"] <- -mypca$var$coord / sqrt(mypca$eig[1, 1])
    }

    # Contribution of original variables to PC1
    mytable[rownames(xcluster), "contribution"] <- mypca$var$contrib / 100

    # Range in pairwise correlations
    if (nrow(xcluster) > 1) {
      mycor <- stats::cor(t(xcluster))
      mycor[upper.tri(mycor, diag = TRUE)] <- NA
      mytable[rownames(xcluster), "min_cor"] <- min(mycor, na.rm = TRUE)
      mytable[rownames(xcluster), "max_cor"] <- max(mycor, na.rm = TRUE)
      mytable[rownames(xcluster), "mean_cor"] <- mean(mycor, na.rm = TRUE)
    }

    # Using PC1-scores to summarise the group
    if (summary == "centroid") {
      if ((sum(mypca$var$coord > 0) / length(mypca$var$coord)) > 0.5) {
        xsummarised[i, ] <- mypca$ind$coord[, 1]
      } else {
        xsummarised[i, ] <- -mypca$ind$coord[, 1]
      }
    }

    # Using the variable most contributing to PC1-scores for summary
    if (summary == "medoid") {
      xsummarised[i, ] <- xcluster[which.max(abs(mypca$var$coord)), ]
      myrownames[i] <- rownames(xcluster)[which.max(abs(mypca$var$coord))]
    }

    if (verbose) {
      utils::setTxtProgressBar(pb, i / k)
    }
  }

  # Re-ordering the description table (multivariate groups first)
  mytable <- mytable[sort.list(mytable$group, decreasing = TRUE), ]

  # Returning the data with summarised features as columns
  colnames(xsummarised) <- colnames(x)
  if (summary == "centroid") {
    rownames(xsummarised) <- paste0("C", 1:k)
  } else {
    rownames(xsummarised) <- myrownames
  }
  xsummarised <- t(xsummarised)
  xsummarised <- scale(xsummarised)

  return(list(xsummarised = xsummarised, description = mytable))
}


# sPLS --------------------------------------------------------------------

#' Calibrate an sPLS model on out-of-sample correlation (for parLapply)
#'
#' Fits an sPLS model with the `param.index`-th (nX, nY) combination of
#' `param` and evaluates it with perf_cor().
#'
#' @param param.index row of `param` to use.
#' @param param matrix with columns `nX` and `nY`.
#' @param Xdata,Ydata predictor and outcome matrices.
#' @param designID subject IDs for multilevel sPLS (optional).
#' @param n_subject number of subjects (multilevel case).
#' @param ncomp number of components.
#' @param nrepeat number of training/test splits.
#' @param seed random seed.
#' @return output of perf_cor().
tune_spls_cor <- function(param.index, param, Xdata, Ydata, designID = NULL, n_subject = NULL, ncomp = 1, nrepeat = 5, seed = 101) {
  keepX <- param[param.index, "nX"]
  keepY <- param[param.index, "nY"]
  if (!is.null(designID)) {
    id.vec <- 1:n_subject
    spls.all <- spls(Xdata, Ydata,
      mode = "regression", keepX = keepX, keepY = keepY, ncomp = ncomp,
      multilevel = designID
    )
  } else {
    id.vec <- NULL
    spls.all <- spls(Xdata, Ydata, mode = "regression", keepX = keepX, keepY = keepY, ncomp = ncomp)
  }

  perf_res <- perf_cor(spls.all, designID = designID, id.vec = id.vec, seed = seed, nrepeat = nrepeat)
  return(perf_res)
}

#' Out-of-sample performance of an sPLS model
#'
#' Repeatedly refits the model on 80% of the subjects (or samples) and
#' computes the correlation between the first X and Y components on the
#' remaining 20%.
#'
#' @param object an sPLS model (mixOmics::spls).
#' @param designID subject IDs for multilevel sPLS (optional).
#' @param id.vec vector of subject IDs to sample from (multilevel case).
#' @param seed random seed.
#' @param nrepeat number of training/test splits.
#' @param keep_value also return the correlation of each split.
#' @return list with the mean and margin of error of the training and test
#'   correlations, VIP, loadings and selection proportions across splits.
perf_cor <- function(object, designID = NULL, id.vec = NULL, seed = 101, nrepeat = 2, keep_value = FALSE) {
  set.seed(seed)
  output <- list()
  seeds_sapply <- sample(1:1000, nrepeat)

  for (i in 1:length(seeds_sapply)) {
    set.seed(seeds_sapply[i]) # set seed for iteration
    if (!is.null(designID)) {
      # Multilevel case: split by subject
      training_ID <- sample(id.vec, length(id.vec) * 0.8)
      test_ID <- id.vec[-which(id.vec %in% training_ID)]

      designID.train <- designID[designID %in% training_ID]
      designID.test <- designID[designID %in% test_ID]
      designID.index.train <- which(designID %in% training_ID)
      designID.index.test <- which(designID %in% test_ID)
      X.train <- object$X[designID.index.train, ]
      X.test <- object$X[designID.index.test, ]
      Y.train <- object$Y[designID.index.train, ]
      Y.test <- object$Y[designID.index.test, ]

      spls.cv <- spls(X.train, Y.train,
        mode = "regression", keepX = object$keepX, keepY = object$keepY, ncomp = object$ncomp,
        multilevel = designID.train
      )
    } else {
      # No multilevel: split by sample
      n <- nrow(object$X)
      row_index <- 1:n
      training_ID <- sample(1:n, n * 0.8)
      test_ID <- row_index[-which(row_index %in% training_ID)]

      X.train <- object$X[training_ID, ]
      X.test <- object$X[test_ID, ]
      Y.train <- object$Y[training_ID, ]
      Y.test <- object$Y[test_ID, ]

      spls.cv <- spls(X.train, Y.train, mode = "regression", keepX = object$keepX, keepY = object$keepY, ncomp = object$ncomp)
    }

    # Extract loadings and scores
    t.cv <- spls.cv$variates$X[, 1]
    u.cv <- spls.cv$variates$Y[, 1]
    a.cv <- spls.cv$loadings$X[, 1]
    b.cv <- spls.cv$loadings$Y[, 1, drop = FALSE]

    # Calculate scores on test data
    t.pred <- X.test %*% a.cv
    u.pred <- Y.test %*% b.cv

    cor_tu_train <- cor(t.cv, u.cv)
    cor_tu_test <- cor(t.pred, u.pred)
    cor_tt_pred <- cor(t.pred, object$variates$X[rownames(t.pred), 1]) # predicted X score and X score with all data
    cor_uu_pred <- cor(u.pred, object$variates$Y[rownames(u.pred), 1]) # predicted Y score and Y score with all data

    output_vec <- c(cor_tu_train, cor_tu_test, cor_tt_pred, cor_uu_pred)
    names(output_vec) <- c("cor_tu_train", "cor_tu_test", "cor_tt_pred", "cor_uu_pred")

    output.temp <- list(cor = output_vec, vip = vip(spls.cv), xloadings = spls.cv$loadings$X, yloadings = spls.cv$loadings$Y)
    output[[i]] <- output.temp
  }

  # Correlations across iterations
  output.cor <- do.call(rbind, lapply(output, function(x) x$cor))

  # Summarise output across iterations
  output.summarised <- list(
    nX = object$keepX[1],
    nY = object$keepY[1],
    cor.tu = cor(object$variates$X[, 1], object$variates$Y[, 1]),
    cor.tu.train.mean = mean(output.cor[, "cor_tu_train"]),
    cor.tu.train.moe = moe(output.cor[, "cor_tu_train"]),
    cor.tu.test.mean = mean(output.cor[, "cor_tu_test"]),
    cor.tu.test.moe = moe(output.cor[, "cor_tu_test"]),
    vip.mean = apply(do.call(rbind, lapply(output, function(x) t(x$vip))), 2, mean),
    vip.moe = apply(do.call(rbind, lapply(output, function(x) t(x$vip))), 2, moe),
    xloadings.mean = apply(do.call(rbind, lapply(output, function(x) t(x$xloadings))), 2, mean),
    yloadings.mean = apply(do.call(rbind, lapply(output, function(x) t(x$yloadings))), 2, mean),
    xloadings.moe = apply(do.call(rbind, lapply(output, function(x) t(x$xloadings))), 2, moe),
    yloadings.moe = apply(do.call(rbind, lapply(output, function(x) t(x$yloadings))), 2, moe),
    x.selection.proportion = apply(do.call(rbind, lapply(output, function(x) t(x$xloadings))), 2, function(x) sum(x != 0) / length(x)),
    y.selection.proportion = apply(do.call(rbind, lapply(output, function(x) t(x$yloadings))), 2, function(x) sum(x != 0) / length(x))
  )

  if (keep_value) {
    output.summarised$cor.train <- output.cor[, "cor_tu_train"]
    output.summarised$cor.test <- output.cor[, "cor_tu_test"]
  }

  return(output.summarised)
}

#' Margin of error of the mean (t-distribution)
#'
#' @param dat numeric vector.
#' @param alpha confidence level.
#' @return half-width of the confidence interval of the mean.
moe <- function(dat, alpha = 0.95) {
  n <- sum(!is.na(dat))
  dat_moe <- qt(1 - (1 - alpha) / 2, n - 1) * sd(dat, na.rm = TRUE) / sqrt(n)
  return(dat_moe)
}

#' Calibrate an sPLS model on cross-validated Q2 and MSEP (for parLapply)
#'
#' @param param.index row of `param` to use.
#' @param param matrix with columns `nX`, `nY` and `comp`.
#' @param Xdata,Ydata predictor and outcome matrices.
#' @param ID subject IDs (multilevel sPLS).
#' @param nrepeat number of repeats of the 5-fold cross-validation.
#' @return named vector: nX, nY, mean Q2 and mean MSEP for each component.
pls.perf.par <- function(param.index, param, Xdata, Ydata, ID, nrepeat = 500) {
  nX <- param[param.index, "nX"]
  nY <- param[param.index, "nY"]
  comp <- param[param.index, "comp"]

  tmp.spls <- spls(Xdata, Ydata,
    ncomp = comp, mode = "regression",
    keepX = nX, keepY = nY, scale = TRUE, multilevel = ID
  )
  tmp.perf <- perf(tmp.spls, validation = "Mfold", folds = 5, nrepeat = nrepeat, criterion = "all", progressBar = FALSE)
  Q2.total.mean.by.comp <- tmp.perf$measures$Q2.total$summary$mean
  MSEP.mean.by.comp <- tmp.perf$measures$MSEP$summary %>%
    group_by(comp) %>%
    summarise(meanMSEP = mean(mean)) %>%
    pull(meanMSEP)
  res <- c(nX, nY, Q2.total.mean.by.comp, MSEP.mean.by.comp)
  names(res) <- c("nX", "nY", paste0("Q2_comp", 1:length(Q2.total.mean.by.comp)), paste0("MSEP_comp", 1:length(MSEP.mean.by.comp)))
  return(res)
}

#' Proportion of variance in the selected Y explained by the X component(s)
#'
#' @param model an sPLS model (mixOmics::spls).
#' @param X,Y predictor and outcome matrices used to fit the model.
#' @return list with the (cumulative) proportion of variance explained, overall
#'   and for each Y variable, by component.
quality.var <- function(model, X, Y) {
  p1 <- predict(model, X)
  matmean <- matrix(colMeans(Y), ncol = dim(Y)[2], nrow = dim(Y)[1], byrow = TRUE)
  sumYT <- colSums((Y - matmean)**2)
  cumP <- rep(0, model$ncomp)
  cumPY <- matrix(0, nrow = model$ncomp, ncol = dim(Y)[2])
  for (k in 1:model$ncomp) {
    selectedY <- which(model$loadings$Y[, k] != 0)
    sumPLS <- colSums((p1$predict[, selectedY, k] - matmean[selectedY])**2)
    cumPY[k, selectedY] <- sumPLS / sumYT[selectedY]
    cumP[k] <- mean(sumPLS / sumYT[selectedY])
  }
  IndP <- diff(c(0, cumP))
  IndPY <- apply(cumPY, 2, FUN = function(x) diff(c(0, x)))
  if (dim(cumPY)[1] > 1) {
    colnames(cumPY) <- colnames(IndPY) <- colnames(Y)
    rownames(cumPY) <- rownames(IndPY) <- paste("comp", 1:model$ncomp, sep = "")
    names(cumP) <- names(IndP) <- paste("comp", 1:model$ncomp, sep = "")
  } else {
    colnames(cumPY) <- names(IndPY) <- colnames(Y)
    names(cumP) <- names(IndP) <- paste("comp", 1:model$ncomp, sep = "")
  }
  return(list(CumPerExplain = cumP, PerExplain = IndP, CumPerExplainY = cumPY, PerExplainY = IndPY))
}


# Validation --------------------------------------------------------------

#' Match metabolic features between two datasets
#'
#' Features are identified as "mass@retention_time". Two features match when
#' their masses differ by less than `ppm` parts per million and their retention
#' times by less than `rt` minutes.
#'
#' @param metab.old,metab.new vectors of feature IDs.
#' @param ppm mass tolerance (ppm).
#' @param rt retention time tolerance (min).
#' @return two-column matrix of matched feature IDs (`old_samemass_v`,
#'   `new_samemass_v`).
compare_metab <- function(metab.old, metab.new, ppm = 7, rt = 0.1) {
  metab.old <- data.frame(OMIC_ID = str_split_fixed(metab.old, ":", 2)[, 1])
  metab.new <- data.frame(OMIC_ID = str_split_fixed(metab.new, ":", 2)[, 1])
  metab.old$mass <- as.numeric(str_split_fixed(metab.old[, "OMIC_ID"], "@", 2)[, 1])
  metab.old$rt <- as.numeric(str_split_fixed(metab.old[, "OMIC_ID"], "@", 2)[, 2])
  metab.new$mass <- as.numeric(str_split_fixed(metab.new[, "OMIC_ID"], "@", 2)[, 1])
  metab.new$rt <- as.numeric(str_split_fixed(metab.new[, "OMIC_ID"], "@", 2)[, 2])

  metab.old <- metab.old[!duplicated(metab.old$OMIC_ID), c("OMIC_ID", "mass", "rt")]
  metab.new <- metab.new[!duplicated(metab.new$OMIC_ID), c("OMIC_ID", "mass", "rt")]

  old_samemass_v <- vector()
  new_samemass_v <- vector()

  ppm <- ppm * 10^(-6)

  for (i in 1:nrow(metab.old)) {
    for (j in 1:nrow(metab.new)) {
      delta <- metab.old$mass[i] * ppm
      if ((metab.old$mass[i] > metab.new$mass[j] - delta & metab.old$mass[i] < metab.new$mass[j] + delta) &
        (metab.old$rt[i] > metab.new$rt[j] - rt & metab.old$rt[i] < metab.new$rt[j] + rt)) {
        old_samemass_v <- c(old_samemass_v, metab.old$OMIC_ID[i])
        new_samemass_v <- c(new_samemass_v, metab.new$OMIC_ID[j])
      }
    }
  }

  return(cbind(old_samemass_v, new_samemass_v))
}
