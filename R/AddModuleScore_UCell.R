#' Calculate module enrichment scores from single-cell data (Seurat v4-only)
#'
#' This v4-only version assumes Seurat::GetAssayData uses the 'slot' argument.
#' @export
AddModuleScore_UCell_v4 <- function(obj, features, maxRank=1500,
        chunk.size=100, BPPARAM=NULL, ncores=1, storeRanks=FALSE,
        w_neg=1, assay=NULL, slot="counts", ties.method="average",
        missing_genes = c("impute","skip"),
        force.gc=FALSE, name="_UCell") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Function 'AddModuleScore_UCell' requires the Seurat package.
            Please install it.", call. = FALSE)
    }
    features <- check_signature_names(features)
    missing_genes <- match.arg(missing_genes)

    if (is.null(assay)) {
        assay <- Seurat::DefaultAssay(obj)
    }

    # If rank matrix was pre-computed in an assay named "UCellRanks", use it
    if ("UCellRanks" %in% Seurat::Assays(obj)) {
        ranks_mat <- tryCatch({
            Seurat::GetAssayData(obj, assay="UCellRanks", slot=slot)
        }, error = function(e) {
            # fallback: try GetAssayData without slot
            tryCatch({
                Seurat::GetAssayData(obj, assay="UCellRanks")
            }, error = function(e2) {
                stop("Could not retrieve 'UCellRanks' assay data: ", conditionMessage(e2))
            })
        })

        meta.list <- rankings2Uscore(
            ranks_mat,
            features=features, chunk.size=chunk.size, w_neg=w_neg,
            missing_genes=missing_genes,
            ncores=ncores, BPPARAM=BPPARAM, force.gc=force.gc, name=name)

    } else {
        # v4-style: try common slots. We'll try the provided slot first, then fallback.
        candidate_slots <- unique(c(slot, "counts", "data", "scale.data"))
        # collect successes
        successful <- list()
        for (s in candidate_slots) {
            dat_try <- tryCatch({
                Seurat::GetAssayData(obj, assay=assay, slot=s)
            }, error = function(e) NULL)
            if (!is.null(dat_try) && (is.matrix(dat_try) || inherits(dat_try, "dgCMatrix") || inherits(dat_try, "Matrix"))) {
                successful[[s]] <- dat_try
            }
        }

        if (length(successful) == 0) {
            stop(sprintf("Cannot find slot '%s' (or any common slots) in assay '%s'.", slot, assay))
        }

        meta.list <- lapply(names(successful)[1], function(x) {
            calculate_Uscore(
                successful[[x]],
                features=features, maxRank=maxRank,
                chunk.size=chunk.size, w_neg=w_neg,
                ncores=ncores, BPPARAM=BPPARAM, ties.method=ties.method,
                missing_genes=missing_genes,
                force.gc=force.gc, storeRanks=storeRanks, name=name)
        })

        meta.list <- unlist(meta.list, recursive = FALSE)

        if (storeRanks==TRUE){
            cells_rankings.merge <- lapply(meta.list,
                function(x) rbind(x[["cells_rankings"]]))
            cells_rankings.merge <- Reduce(cbind, cells_rankings.merge)
            obj[["UCellRanks"]] <- Seurat::CreateAssayObject(cells_rankings.merge)
        }
    }

    meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_U"]]))
    meta.merge <- Reduce(rbind, meta.merge)
    obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))
    return(obj)
}
