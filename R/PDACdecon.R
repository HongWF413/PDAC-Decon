#' Transform protein name to gene symbol
#'
#' @param x character vector, protein id used to be mapped to gene
#' @param org.db org.db object
#' @param match_table a data.frame which must have two columns: FROM, TO
#' @param org.db.from keytypes which annotationdbi used
#' @param org.db.to keytypes which annotationdbi used
#'
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keytypes
#' @export
#' @examples
#' # protein2gene()
#'

protein2gene <- function(
    x,
    org.db,
    match_table,
    org.db.from = "ALIAS",
    org.db.to = "SYMBOL"
){ # browser()

  stopifnot(is.vector(x) && !is.list(x))

  method_var <- "none"

  if(missing(org.db)){
    if(missing(match_table)){
      stop("org.db or match_table should be provided")
    }else{
      method_var <- "match_table"
    }
  }else{
    if(missing(match_table)){
      inherits(org.db, "OrgDb")
      method_var <- "org.db"
    }else{
      warning("org.db and match_table provided, only used match_table", immediate. = TRUE)
      method_var <- "match_table"
    }
  }

  # 正式转换
  if(method_var == "org.db"){

    orgdb_keytypes <-  keytypes(org.db)

    stopifnot(org.db.from %in% orgdb_keytypes)
    stopifnot(org.db.to %in% orgdb_keytypes)

    res <-
      select(
        org.db,
        keys = x,
        columns = org.db.to,
        keytype = org.db.from
      )

    protein_matched_count <- sum(unique(res[[org.db.from]]) %in% x)
    if(protein_matched_count < 1){
      stop("None of protein can be matched using org.db")
    }else if(protein_matched_count/length(x) < 0.5){
      warning("Matched protein less than 50%")
    }

    return(res)
  }

  # method_var: match_table
  stopifnot(all(c("FROM", "TO") %in% colnames(match_table)))
  if(length(x) > nrow(match_table)){
    warning("row length of match_table less than x, please check you match_table", immediate. = TRUE)
  }


  x_matched <- match_table$TO[match(x, match_table$FROM)]

  res <- data.frame(FROM = x, TO = x_matched)

  protein_matched_count <- sum(!is.na(x_matched))
  if(protein_matched_count < 1){
    stop("None of protein can be matched using org.db")
  }else if(protein_matched_count/length(x) < 0.5){
    warning("Matched protein less than 50%")
  }

  return(res)

}

if(FALSE){
  celldive_mean_exprs$CellType %>% protein2gene(org.Hs.eg.db)
  celldive_mean_exprs$CellType %>% protein2gene(match_table = ll %>% `colnames<-`(c("FROM", "TO")))
}




#' Create Profile Matrix for SpatialDecon from CellDive Data.
#'
#' This function generate a Profile Matrix for SpatialDecon from CellDive Data.
#'
#' @param x celldive average expression dataset, a data.frame or a matrix
#' @param bulk_counts bulk RNA expression dataset, a data.frame or a matrix
#' @param celltypes a data.frame which must have two columns: gene, cell
#' @param protein_to_gene logical value, if TURE rownames of celldive will be used as protein, and will be mapped to gene id
#' @param correlation most similar gene cutoff
#' @param weight_factor weight value of a gene when it was a celltype marker
#' @param verbose if TRUE, then more message will be print
#' @param org.db org.db object, only use when protein_to_gene is TURE
#' @param match_table a data.frame which must have two columns: FROM, TO, only use when protein_to_gene is TURE
#' @param org.db.from keytypes which annotationdbi used, only use when protein_to_gene is TURE
#' @param org.db.to keytypes which annotationdbi used, only use when protein_to_gene is TURE
#' @importFrom stats cor
#' @export
#' @examples
#' # PDAC_Decon()
#'

PDAC_Decon <- function(
    x, # celldive average expr, data.frame matrix
    bulk_counts, # data.frame matrix
    celltypes, # gene, protein, data.frame matrix
    protein_to_gene = TRUE,
    correlation = 0.4,
    weight_factor = 0.2,
    verbose = TRUE,

    org.db,
    match_table,
    org.db.from = "ALIAS",
    org.db.to = "SYMBOL"
){
  # browser()

  # x: data.frame
  # rownames: protein
  # column:cell
  stopifnot(is.data.frame(x) || is.matrix(x))
  x <- as.data.frame(x)

  is_all_numeric <- lapply(x, function(x) is.numeric(x) )
  if(!all(unlist(is_all_numeric))){
    stop("All values in celldive average expression dataset should be numeric.")
  }

  # celldive potein -> gene
  if(protein_to_gene){
    if(verbose) message("Protein name will be mapped to gene Symbol.")

    gene_table <- protein2gene(
      rownames(x),
      org.db      = org.db,
      match_table = match_table,
      org.db.from = org.db.from,
      org.db.to   = org.db.to
    )
    gene_table <- gene_table[!is.na(gene_table[[2]]), ,drop=FALSE]

    # protein expr -> gene expr
    x <- x[gene_table[[1]], ]
    rownames(x) <- gene_table[[2]]

    # TODO: gene_table[[2]]可能有重复基因，导致x的行名有问题
  }

  # cor table: celldive gene -> most similar bulk gene
  stopifnot(is.data.frame(bulk_counts) || is.matrix(bulk_counts))
  bulk_counts <- as.data.frame(bulk_counts)

  is_all_numeric2 <- lapply(bulk_counts, function(x) is.numeric(x) )
  if(!all(unlist(is_all_numeric2))){
    stop("All values in bulk expression dataset should be numeric.")
  }

  # gene -> most similar gene
  cd_gene <- rownames(x)
  bulk_gene <- rownames(bulk_counts)
  cd_gene_exist_in_bulk_counts <- intersect(cd_gene, bulk_gene)

  if(length(cd_gene_exist_in_bulk_counts) < 1){
    if(verbose){
      if(protein_to_gene){
        cat("############ Verbose message ###########")
        cat("Protein had been mapped to gene: \n")
        print(gene_table)
      }else{
        cat("############ Verbose message ###########")
        cat("Protein not in bulk counts: \n")
        print(cd_gene)
      }
    }
    stop("None of celldive gene existed in bulk counts.")
  }


  bulk_counts_sub <- bulk_counts[cd_gene_exist_in_bulk_counts, ,drop=FALSE]

  cor_table <- cor(t(bulk_counts_sub), t(bulk_counts))
  cor_ls_filter <- apply(cor_table, 1, function(x) rownames(bulk_counts)[x > correlation], simplify = FALSE)

  most_similar_gene <- unique(unlist(cor_ls_filter, use.names = FALSE))

  # 先对ratio2做必要的处理
  if(!missing(celltypes)){
    stopifnot(is.data.frame(celltypes) || is.matrix(celltypes))
    if(!all(colnames(celltypes) %in% c("cell", "gene"))){
      stop(
        "celltypes should have two columns: cell and gene"
      )
    }

    cell_in_celldive <- unique(celltypes[['cell']]) %in% colnames(x)
    if(sum(cell_in_celldive) < 1){
      stop("None of cell of celltypes found in celldive dataset.")
    }else if(sum(cell_in_celldive) < length(unique(celltypes[['cell']]))){
      warning("Not all cell of celltypes found in celldive dataset.")
    }

    gene_in_celldive <- unique(celltypes[['gene']]) %in% rownames(x)
    if(sum(gene_in_celldive) < 1){
      stop("None of gene of celltypes found in celldive dataset.")
    }else if(sum(gene_in_celldive) < length(unique(celltypes[['cell']]))){
      warning("Not all gene of celltypes found in celldive dataset.")
    }

    # celltypes <- celltypes[gene_in_celldive, cell_in_celldive, drop=FALSE]
    celltypes <- celltypes[celltypes$gene %in% rownames(x) & celltypes$cell %in% colnames(x), , drop = FALSE]

  }else{
    #TODO: 很可能导致后续的SpatialDecon报错
  }


  all_celldive_celltypes <- colnames(x)

  res_ls <- list()
  for(id in seq_along(cor_ls_filter)){
    gene1 <- names(cor_ls_filter) [id]
    gene2 <- cor_ls_filter [[id]]

    # WTA -> mean per mean
    bulk_expr_sub <- bulk_counts[gene2, ,drop=FALSE]
    bulk_mean_expr_sub <- apply(bulk_expr_sub, 1, mean)

    # ratio1
    celldive_sub <- unlist(x[gene1, ])
    celldive_sub_ratio <- celldive_sub/sum(celldive_sub)

    celldive_expand_sub <- matrix(0, nrow = length(gene2), ncol = length(all_celldive_celltypes))
    celldive_expand_sub <- sweep(celldive_expand_sub, 2, celldive_sub_ratio, FUN = "+")

    ratio1 <- celldive_expand_sub

    # ratio2

    # fix ratio2 according to celltypes
    # if(id == 3) browser()
    if(
      !missing(celltypes) &&
      any(celltypes$gene %in% gene2)
    ){
      ratio2 <- matrix(
        weight_factor,
        nrow = length(gene2),
        ncol = length(all_celldive_celltypes),
        dimnames = list(
          gene2,
          all_celldive_celltypes
        )
      )
      # subset
      gene_idx <- intersect(gene2, celltypes$gene)
      cell_idx <- intersect(all_celldive_celltypes, celltypes$cell)

      ratio2[gene_idx, cell_idx] <- 1- weight_factor

    }else{
      ratio2 <- matrix(
        1,
        nrow = length(gene2),
        ncol = length(all_celldive_celltypes),
        dimnames = list(
          gene2,
          all_celldive_celltypes
        )
      )
    }

    # result
    res <- bulk_mean_expr_sub * ratio1 * ratio2
    res_df_sub <- as.data.frame(res)
    res_df_sub$gene <- rownames(res)

    # res_ls[[g1]] <- res_df_sub
    res_ls[[id]] <- res_df_sub
  }

  res_df <- do.call(rbind, res_ls)

  # average cell expression grouped by gene column(the last column is gene)
  res_df_remove_gene <- res_df[colnames(res_df) != "gene"]
  res_final <- sapply(res_df_remove_gene, function(x) tapply(x, res_df$gene, mean))

  res_final
}

if(FALSE){
  PDAC_Decon(
    x = celldive_mean_exprs,
    protein_to_gene = F,
    bulk_counts = counts_ana
  )

  # celltypes
  celltypes_test <-
    data.frame(
      cell = c("B.cell","B.cell.PD.L1","B.cell.PD.L1","Dendritic.Cells","Dendritic.Cells","Dendritic.Cells","Dendritic.Cells.PD.1","Dendritic.Cells.PD.1","Dendritic.Cells.PD.1","Dendritic.Cells.PD.1",
               "Fibroblast","Fibroblast","Fibroblast","Fibroblast","Fibroblast","M1.Macrophage","M1.Macrophage","M1.Macrophage","M1.Macrophage","M1.Macrophage","M1.Macrophage"
      ),
      gene = c("CD3G","CD14","MKI67","CD8A","CD3G","FOXP3","CD33","ITGAM",
               "SNCA","SMN1","MKI67","SPATA2","CD33","KRT1","MS4A1","CD274","CD14","CD68","CD8A","CD163","CD33"
      )
    )

  PDAC_Decon(
    x = celldive_mean_exprs,
    protein_to_gene = F,
    celltypes = celltypes_test,
    bulk_counts = counts_ana
  )
}
