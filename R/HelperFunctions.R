#' Calculate Mann Whitney U from a vector of ranks
u_stat <- function(ranks_matrix, gene_idx, maxRank=1500, sparse=FALSE){
  
    len_sig <- length(gene_idx)
    ncells <- ncol(ranks_matrix)-1 
    
    if (len_sig <= 0) {
      return(rep(0, ncells))
    }
    present_idx <- gene_idx[gene_idx>0]
    missing_idx <- gene_idx[gene_idx<0]
  
    #Initialize rank_sum with missing genes
    rank_sum <- rep(length(missing_idx) * maxRank, ncells)
  
    #Subset ranks_matrix by rows (genes) using numeric indices
    if (length(present_idx) > 0) {
      rank_sub <- ranks_matrix[present_idx, -1, drop = FALSE]
      if (sparse==TRUE){
        rank_sub[rank_sub==0] <- maxRank
      }
      insig <- rank_sub >= maxRank
      rank_sub[insig] <- maxRank
    
      rank_sum = rank_sum + apply(rank_sub, 2, sum)
    }
  
    rank_sum_min <- len_sig*(len_sig + 1)/2
    ucell_score <- 1 - (rank_sum-rank_sum_min)/(len_sig*maxRank - rank_sum_min)
    return(ucell_score)
}

#' Calculate U scores for a list of signatures, given a rank matrix
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1500,
    sparse=FALSE, w_neg=1, missing_genes="impute") {
    
    cell_names <- colnames(ranks_matrix)[-1]
    dim <- length(cell_names)
    all_genes <- ranks_matrix[["rn"]]
    
    u_matrix <- vapply(sig_list, FUN.VALUE = numeric(dim), FUN=function(sig) {
        sig_neg <- grep('-$', unlist(sig), perl=TRUE, value=TRUE)
        sig_pos <- setdiff(unlist(sig), sig_neg)
        
        #Positive gene set
        sig_pos <- gsub('\\+$','',sig_pos,perl=TRUE)
        pos_idx <- get_gene_idx(all_genes, sig_pos,
                                 missing_genes = missing_genes)
        
        u_p <- u_stat(ranks_matrix, pos_idx, maxRank, sparse=sparse)

        #Negative gene set
        sig_neg <- gsub('-$','',sig_neg,perl=TRUE)
        neg_idx <- get_gene_idx(all_genes, sig_neg,
                                 missing_genes = missing_genes)
        
        u_n <- u_stat(ranks_matrix, neg_idx, maxRank, sparse=sparse)

        diff <- u_p - w_neg*u_n   #Subtract negative sets, if any
        diff[diff<0] <- 0  #clip negative values
        return(diff)
    })
    if (is.vector(u_matrix)) {  # Case of ncells=1
      u_matrix <- t(as.matrix(u_matrix))
    }
    
    rownames(u_matrix) <- cell_names
    return (u_matrix)
}

#' Calculate rankings and scores for query data and given signature set
calculate_Uscore <- function(
        matrix, features,  maxRank=1500, chunk.size=100,
        BPPARAM = NULL, ncores=1, w_neg=1, ties.method="average",
        missing_genes = c("impute","skip"),
        storeRanks=FALSE, force.gc=FALSE, name="_UCell"){
    
    #Make sure we have a sparse matrix
    if (!methods::is(matrix, "dgCMatrix")) {
        matrix <- Matrix::Matrix(as.matrix(matrix),sparse = TRUE)
    }
    missing_genes <- match.arg(missing_genes)
    
    #Do not evaluate more genes than there are
    if (!is.numeric(maxRank)) {
        stop("Rank cutoff (maxRank) must be a number")
    }
    if (maxRank > nrow(matrix)) {
        maxRank <- nrow(matrix)
    }
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (w_neg<0) {stop("Weight on negative signatures (w_neg) must be >=0")}
    
    #Signatures cannot be larger than maxRank parameter
    sign.lgt <- lapply(features, length)
    if (any(sign.lgt > maxRank)) {
        stop("One or more signatures contain more genes than maxRank parameter.
            Increase maxRank parameter or make shorter signatures")
    }
    
    #Split into manageable chunks
    split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)
    
    #Either take a BPPARAM object, or make one on the spot using 'ncores'
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
    }
    meta.list <- BiocParallel::bplapply(
        X = split.data, 
        BPPARAM =  BPPARAM,
        FUN = function(x) {
            cells_rankings <- data_to_ranks_data_table(x,
                ties.method=ties.method)
            cells_U <- u_stat_signature_list(features, cells_rankings, 
                maxRank=maxRank, sparse=FALSE, missing_genes=missing_genes,
                w_neg=w_neg)
            colnames(cells_U) <- paste0(colnames(cells_U),name)
            if (storeRanks==TRUE){
                gene.names <- as.character(as.matrix(cells_rankings[,1]))
                #make sparse (rank=0 means rank>=maxRank)
                cells_rankings[cells_rankings>=maxRank] <- 0
                ranks.sparse <- Matrix::Matrix(as.matrix(
                    cells_rankings[,-1]),sparse = TRUE)
                dimnames(ranks.sparse)[[1]] <- gene.names
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_rankings=ranks.sparse, cells_U=cells_U))
            } else {
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_U=cells_U))
            }
            
        })
    return(meta.list)
}

#' Get signature scores from pre-computed rank matrix
rankings2Uscore <- function(ranks_matrix, features, chunk.size=100, w_neg=1,
                            BPPARAM = NULL,ncores=1, force.gc=FALSE,
                            missing_genes = c("impute","skip"),
                            name="_UCell") {
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (!is.numeric(w_neg) | w_neg<0) {
        stop("Weight on negative signatures (w_neg) must be >=0")}
    
    maxRank <- max(ranks_matrix)+1
    split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
    rm(ranks_matrix)
    
    #Either take a BPPARAM object, or make one on the spot using 'ncores'
    if (is.null(BPPARAM)) {
        BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
    }
    meta.list <- BiocParallel::bplapply(
        X = split.data, 
        BPPARAM =  BPPARAM,
        FUN = function(x) {
            
            dense <- as.matrix(x)
            dense <- as.data.table(dense, keep.rownames=TRUE)
            setkey(dense, "rn", physical=FALSE)
            
            cells_U <- u_stat_signature_list(features, dense,
                maxRank=maxRank, sparse=TRUE, 
                missing_genes = missing_genes, w_neg=w_neg)
            colnames(cells_U) <- paste0(colnames(cells_U),name)
            
            if (force.gc) {
                dense <- NULL
                gc()
            }
            return(list(cells_U=cells_U))
        }
    )
    return(meta.list)
}

#' Get indices for genes in data matrix
get_gene_idx <- function(all_genes, signature, missing_genes="impute") {
  
  idx <- match(signature, all_genes)  # get numeric indices
  
  if (missing_genes == "skip") {
    idx <- idx[!is.na(idx)]
  } else if (missing_genes == "impute") {
    # Replace missing genes with -1; we'll handle them later in scoring
    idx[is.na(idx)] <- -1
  }
  return(idx)
}

#' Check signature names and add standard names is missing
check_signature_names <- function(features) {
    defaultSigName <- paste0(rep("signature_",length(features)),
        seq_along(features))
    if(is.null(names(features))){
        names(features) <- defaultSigName
    } else {
        invalidNames <- names(features) == "" | duplicated(names(features))
        names(features)[invalidNames] <- defaultSigName[invalidNames]
    }
    return(features)
}

#' Calculate per-cell feature rankings
data_to_ranks_data_table <- function(data, ties.method="average") {
    dt <- as.data.table(as.matrix(data))
    rnaDT.ranks.dt <- dt[, lapply(.SD, function(x)
        frankv(x,ties.method=ties.method,order=c(-1L)))]
    rnaDT.ranks.rownames <- rownames(data)
    rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
    setkey(rnaDT.ranks.dt.rn, "rn", physical = FALSE)
    return(rnaDT.ranks.dt.rn)
}

#' Split data matrix into smaller sub-matrices ('chunks')
split_data.matrix <- function(matrix, chunk.size=100) {
    ncols <- dim(matrix)[2]
    nchunks <- (ncols-1) %/% chunk.size + 1
    
    split.data <- list()
    min <- 1
    for (i in seq_len(nchunks)) {
        if (i == nchunks-1) {  #make last two chunks of equal size
            left <- ncols-(i-1)*chunk.size
            max <- min+round(left/2)-1
        } else {
            max <- min(i*chunk.size, ncols)
        }
        split.data[[i]] <- matrix[,min:max,drop=FALSE]
        min <- max+1    #for next chunk
    }
    return(split.data)
}

#' Smoothing scores by KNN
knn_smooth_scores <- function(
    matrix=NULL,
    nn=NULL,
    decay=0.1,   #decay must be bound between 0 and 1
    up.only=FALSE #scores can only increase
) {
  
  sig.cols <- colnames(matrix)
  
  w.df <- vapply(sig.cols, FUN.VALUE=numeric(nrow(matrix)), FUN=function(s) {
    ss.scores <- matrix[,s]
    weighted.scores <- vapply(X = seq_len(nrow(nn$index)),
                              FUN.VALUE = numeric(1),
                              FUN = function(x) {
                                r <- nn$index[x,]
                                r <- c(x,r)
                                i <- seq(0, length(r)-1)
                                w <- (1-decay)**i
                                sum(w * ss.scores[r])/sum(w)
                              })
    if (up.only) {
      pmax(weighted.scores, ss.scores)
    } else {
      weighted.scores
    }
  })
  rownames(w.df) <- rownames(matrix)
  as.data.frame(w.df)
}  

#' @rdname SmoothKNN
#' @method SmoothKNN Seurat
#' @export
SmoothKNN.Seurat <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="pca",
    k=10,
    decay=0.1,
    up.only=FALSE,
    BNPARAM=AnnoyParam(),
    BPPARAM=SerialParam(),
    suffix="_kNN",
    assay=NULL,
    slot="data",
    sce.expname=NULL,
    sce.assay=NULL
) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Function 'SmoothKNN_UCell' requires the Seurat package.
            Please install it.", call. = FALSE)
  } 
  if (!reduction %in% Seurat::Reductions(obj)) {
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  if (is.null(signature.names)) {
    stop("Please provide the metadata column names that you want to smooth")
  }
  
  if (is.null(assay)) {  # Work on metadata
    found <- intersect(signature.names, colnames(obj[[]]))
    notfound <- setdiff(signature.names, found)
    
    if (length(found)==0) {
      stop("Could not find any of the given signatures in this object")
    }
    if (length(notfound)>0) {
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following signature were found in metadata:\n* %s",nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- obj[[found]]
  } else {  # Work directly on features
    exp <- Seurat::GetAssayData(obj, slot=slot, assay=assay)
    feats <- rownames(exp)
    found <- intersect(signature.names, feats)
    notfound <- setdiff(signature.names, found)
    
    if (length(found)==0) {
      stop("Could not find any of the given features in this object")
    }
    if (length(notfound)>0) {
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following features were not found in assay %s:\n* %s",
                      assay, nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- t(exp[found, , drop=FALSE])
  }
  ncells <- ncol(obj)
  
  if (decay<0 | decay>1) {
    stop("decay parameter must be a number between 0 and 1")
  }
    
  if (k<=0) {  #this behavior disables kNN smoothing
    k=1
    decay=1
  }
  
  if (ncells <= k) {
    k <- ncells-1
    warning("'k' capped at the number of observations minus 1")
  }
  
  if (ncells>1) {
    # Find kNNs
    space <- Seurat::Embeddings(obj, reduction=reduction)
    nn <- findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    
    # Do smoothing
    smooth.df <- knn_smooth_scores(matrix=m, nn=nn,
                                   decay=decay, up.only=up.only)  
  } else {
    smooth.df <- m
  }
  
  if (is.null(assay)) {  #metadata
    colnames(smooth.df) <- paste0(colnames(smooth.df), suffix)
    obj <- Seurat::AddMetaData(obj, metadata = smooth.df)
  } else {  #new assay
    nas <- paste0(assay, suffix)
    obj[[nas]] <- Seurat::CreateAssayObject(data=t(smooth.df))
  }
  return(obj)
}
