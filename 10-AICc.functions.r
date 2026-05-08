get.params <- function (model) {
    lambdas <- model[1:(length(model) - 4)]
    countNonZeroParams <- function(x) {
        if (strsplit(x, split = ", ")[[1]][2] != "0.0") 
            1
    }
    no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
    return(no.params)
}

# adapted
calc.aicc <- function (nparam, occ, predictive.maps) {
    AIC.valid <- nparam < length(occ)
		vals <- occ
        probsum <- sum(predictive.maps)
 
        LL <- colSums(log(t(t(vals)/probsum)), na.rm = TRUE)
        AICc <- (2 * nparam - 2 * LL) + (2 * (nparam) * (nparam + 
            1)/(length(occ) - nparam - 1))
        AICc[AIC.valid == FALSE] <- NA
        AICc[is.infinite(AICc)] <- NA
        if (sum(is.na(AICc)) == length(AICc)) {
            warning("AICc not valid... returning NA's.")
            res <- data.frame(cbind(AICc, delta.AICc = NA, w.AIC = NA, 
                parameters = nparam))
        }
        else {
            delta.AICc <- (AICc - min(AICc, na.rm = TRUE))
            w.AIC <- (exp(-0.5 * delta.AICc))/(sum(exp(-0.5 * 
                delta.AICc), na.rm = TRUE))
            res <- data.frame(AICc, delta.AICc, w.AIC, parameters = nparam)
            rownames(res) <- NULL

    }
    rownames(res) <- NULL
    return(res)
}





