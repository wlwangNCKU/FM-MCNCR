# classError function in package: mclust
classError = function(classification, class)
{
    q <- function(map, len, x) {
        x <- as.character(x)
        map <- lapply(map, as.character)
        y <- sapply(map, function(x) x[1])
        best <- y != x
        if (all(len) == 1)
            return(best)
        errmin <- sum(as.numeric(best))
        z <- sapply(map, function(x) x[length(x)])
        mask <- len != 1
        counter <- rep(0, length(len))
        k <- sum(as.numeric(mask))
        j <- 0
        while (y != z) {
            i <- k - j
            m <- mask[i]
            counter[m] <- (counter[m]%%len[m]) + 1
            y[x == names(map)[m]] <- map[[m]][counter[m]]
            temp <- y != x
            err <- sum(as.numeric(temp))
            if (err < errmin) {
                errmin <- err
                best <- temp
            }
            j <- (j + 1)%%k
        }
        best
    }
    if (any(isNA <- is.na(classification))) {
        classification <- as.character(classification)
        nachar <- paste(unique(classification[!isNA]), collapse = "")
        classification[isNA] <- nachar
    }
    MAP <- mapClass(classification, class)
    len <- sapply(MAP[[1]], length)
    if (all(len) == 1) {
        CtoT <- unlist(MAP[[1]])
        I <- match(as.character(classification), names(CtoT),
            nomatch = 0)
        one <- CtoT[I] != class
    }
    else {
        one <- q(MAP[[1]], len, class)
    }
    len <- sapply(MAP[[2]], length)
    if (all(len) == 1) {
        TtoC <- unlist(MAP[[2]])
        I <- match(as.character(class), names(TtoC), nomatch = 0)
        two <- TtoC[I] != classification
    }
    else {
        two <- q(MAP[[2]], len, classification)
    }
    err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
        as.vector(one)
    else as.vector(two)
    bad <- seq(along = classification)[err]
    list(misclassified = bad, errorRate = length(bad)/length(class))
}

mapClass = function(a, b)
{
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
        warning("unequal lengths")
        return(x)
    }
    if (is.factor(a) & is.factor(b) & nlevels(a) == nlevels(b)) {
        aTOb <- as.list(levels(b))
        names(aTOb) <- levels(a)
        bTOa <- as.list(levels(a))
        names(bTOa) <- levels(b)
        out <- list(aTOb = aTOb, bTOa = bTOa)
        return(out)
    }
    if (is.character(a) & is.character(b) & length(unique(a)) ==
        length(unique(b))) {
        aTOb <- as.list(unique(b))
        names(aTOb) <- unique(a)
        bTOa <- as.list(unique(a))
        names(bTOa) <- unique(b)
        out <- list(aTOb = aTOb, bTOa = bTOa)
        return(out)
    }
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
        I <- match(Max[i], Tab[i, ], nomatch = 0)
        aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b))
        aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
        J <- match(Max[j], Tab[, j])
        bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a))
        bTOa <- lapply(bTOa, as.numeric)
    out <- list(aTOb = aTOb, bTOa = bTOa)
    return(out)
}

# adjustedRandIndex function in package: mclust
adjustedRandIndex = function(x, y)
{
    x <- as.vector(x)
    y <- as.vector(y)
    if (length(x) != length(y))
        stop("arguments must be vectors of the same length")
    tab <- table(x, y)
    if (all(dim(tab) == c(1, 1)))
        return(1)
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c
    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    return(ARI)
}
