#' @title Summarize
#' @description Summarize region.
#' @return A list.
summarize.region.parameters <- function (region.genes, model, X, Y, grouping.th = 0.9, rm.na = TRUE) {

  # Take average of the Zs and Ws over the overlapping models
  # Useful for interpretation.
  # region.genes <- top.regs[[i]]; X <- ge; Y <- cn; grouping.th = 0.8
  
  zs <- array(NA, dim = c(ncol(X$data), length(region.genes)))
  wxs <- wys <- array(NA, dim = c(length(region.genes), length(region.genes)))
  colnames(zs) <- colnames(wxs) <- colnames(wys) <- rownames(wxs) <- rownames(wys) <- region.genes
  rownames(zs) <- colnames(X$data)

  chr <- na.omit(as.numeric(as.character(unique(X$info[region.genes,"chr"]))))
  chr.models <- model@chromosomeModels[[chr]]@models
  # FIXME: provide function that automatically fetches model names!
  model.names <- sapply(1:length(chr.models), function (k) {chr.models[[k]]@geneName})
  # Pick genes from this region that also have models
  gs <- intersect(model.names, region.genes)
  # Go through these models
  for (g in gs) {

    m <- findModel(model, g)

    z <- z.effects(m, X, Y)
    zs[rownames(z), g] <- z

    # FIXME: sign switch not taken into account in m@W!
    #wx <- m@W$X[rownames(m@W$X) %in% region.genes,]
    #wy <- m@W$Y[rownames(m@W$Y) %in% region.genes,]
    we <- W.effects(m, X, Y)
    wx <- we$X
    wy <- we$Y
    comx <- intersect(names(wx), rownames(wxs))
    comy <- intersect(names(wy), rownames(wys))
    wxs[comx, g] <- wx[comx]
    wys[comy, g] <- wy[comy]
    
  }

  
  wxs <- wxs[,!apply(wxs, 2, function(x){all(is.na(x))})]
  wys <- wys[,!apply(wys, 2, function(x){all(is.na(x))})]
  zs <- zs[,!apply(zs, 2, function(x){all(is.na(x))})]  
    
  
  # Now it is possible that the region contains multiple different alteration regions
  # go through zs and detect groups of neighboring and highly correlated models
  cors <- (cor(zs) > grouping.th)
  k <- 1
  groups <- list()
  while (length(cors)>0) {
    inds <- which(cors[1, ])
    groups[[k]] <- names(inds)
    k <- k+1
    cors <- cors[-inds, -inds, drop = FALSE]
  }
  
  # take means within each group
  summaries <- list()
  for (k in 1:length(groups)) {
    gs <- groups[[k]]
    wx.tmp <- rowMeans(matrix(wxs[,gs], nrow(wxs)), na.rm = TRUE)
    names(wx.tmp) <- rownames(wxs)
    wy.tmp <- rowMeans(matrix(wys[,gs], nrow(wys)), na.rm = TRUE)
    names(wy.tmp) <- rownames(wys)

    if (rm.na) {
      wx.tmp <- wx.tmp[!is.na(wx.tmp)]
      wy.tmp <- wy.tmp[!is.na(wy.tmp)]
    }
    
    W <- list(X = wx.tmp, Y = wy.tmp)
    Z <- rowMeans(matrix(zs[,gs], nrow(zs)), na.rm = TRUE)
    names(Z) <- rownames(zs)
    summaries[[k]] <- list(z = Z, W = W )
  }

  summaries
  
}


