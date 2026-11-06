# =============================================================================
# object_export.R — slim Seurat .rds and AnnData .h5ad for the customer export
# Sourced by 08f_customer_export.R. No new dependencies: the h5ad is written
# with hdf5r (already in the env) following the AnnData on-disk spec 0.1.0/0.2.0.
# =============================================================================
OBJ_META <- c("orig.ident", "sample", "nCount_RNA", "nFeature_RNA", "percent.mt",
              "seurat_clusters", "cell_type", "cell_type_ref", "cell_type_denovo")

write_h5ad <- function(path, counts, obs, umap = NULL) {
  stopifnot(requireNamespace("hdf5r", quietly = TRUE), inherits(counts, "dgCMatrix"))
  # UTF-8 variable-length strings: with hdf5r's default ASCII charset h5py hands anndata bytes (b'A')
  str_t <- hdf5r::H5T_STRING$new(size = Inf)$set_cset(hdf5r::h5const$H5T_CSET_UTF8)
  # scalar string / bool attributes: AnnData readers expect a scalar dataspace,
  # while hdf5r's `h5attr<-` writes a length-1 array
  sattr <- function(h, nm, v) {
    a <- h$create_attr(nm, dtype = if (is.character(v)) str_t else hdf5r::h5types$H5T_NATIVE_HBOOL, space = hdf5r::H5S$new("scalar"))
    a$write(v); a$close()
  }
  enc <- function(h, type, ver) { sattr(h, "encoding-type", type); sattr(h, "encoding-version", ver) }
  strds <- function(g, nm, v) g$create_dataset(nm, robj = as.character(v), dtype = str_t)
  if (file.exists(path)) file.remove(path)
  f <- hdf5r::H5File$new(path, mode = "w"); on.exit(f$close_all(), add = TRUE)
  enc(f, "anndata", "0.1.0")
  # X: cells x genes CSR. A genes x cells dgCMatrix (CSC) has exactly the arrays of
  # the cells x genes CSR (p = indptr, i = indices, x = data), so no transpose.
  X <- f$create_group("X"); enc(X, "csr_matrix", "0.1.0")
  hdf5r::h5attr(X, "shape") <- c(ncol(counts), nrow(counts))
  xv <- counts@x; if (all(xv == round(xv))) xv <- as.integer(xv)
  X$create_dataset("data", robj = xv)
  X$create_dataset("indices", robj = as.integer(counts@i))
  X$create_dataset("indptr", robj = as.integer(counts@p))
  df <- function(name, index, cols) {
    g <- f$create_group(name); enc(g, "dataframe", "0.2.0")
    sattr(g, "_index", "_index")
    if (length(cols)) hdf5r::h5attr(g, "column-order") <- names(cols) else
      g$create_attr("column-order", dtype = str_t,
                    space = hdf5r::H5S$new("simple", dims = 0, maxdims = 0))$close()
    strds(g, "_index", index)
    for (nm in names(cols)) {
      v <- cols[[nm]]
      if (is.numeric(v) || is.logical(v)) {
        d <- g$create_dataset(nm, robj = if (is.logical(v)) as.integer(v) else v); enc(d, "array", "0.2.0")
      } else {
        v <- as.character(v); lv <- sort(unique(v[!is.na(v)]))
        cg <- g$create_group(nm); enc(cg, "categorical", "0.2.0"); sattr(cg, "ordered", FALSE)
        strds(cg, "categories", lv)
        # AnnData marks a missing value with code -1; widen the integer type past 127 categories
        codes <- match(v, lv) - 1L; codes[is.na(codes)] <- -1L
        cg$create_dataset("codes", robj = as.integer(codes),
                          dtype = if (length(lv) <= 127L) hdf5r::h5types$H5T_NATIVE_INT8 else hdf5r::h5types$H5T_NATIVE_INT32)
      }
    }
  }
  df("obs", colnames(counts), as.list(obs))
  df("var", rownames(counts), list())
  om <- f$create_group("obsm"); enc(om, "dict", "0.1.0")
  if (!is.null(umap)) {
    # hdf5r stores an R matrix with reversed dims, so t() yields (n_obs, 2) in h5py
    d <- om$create_dataset("X_umap", robj = t(unname(as.matrix(umap)))); enc(d, "array", "0.2.0")
  }
  for (g in c("obsp", "varm", "varp", "layers", "uns")) enc(f$create_group(g), "dict", "0.1.0")
  invisible(path)
}

export_objects <- function(obj, md, lab, exp_dir, run, msg) {
  SeuratObject::DefaultAssay(obj) <- "RNA"
  rd <- intersect("umap", SeuratObject::Reductions(obj))
  s  <- Seurat::DietSeurat(obj, assays = "RNA", layers = "counts", dimreducs = rd)
  keep <- intersect(OBJ_META, colnames(s@meta.data))
  s@meta.data <- s@meta.data[, keep, drop = FALSE]
  f <- file.path(exp_dir, paste0(run, "_seurat_slim.rds")); t0 <- Sys.time()
  saveRDS(s, f, compress = "xz")
  msg("%s  %.1f MB  (%.0f s, xz)", basename(f), file.size(f)/1e6,
      as.numeric(Sys.time() - t0, units = "secs"))
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    msg("hdf5r not installed - h5ad skipped"); return(invisible())
  }
  cm <- as(SeuratObject::LayerData(s, layer = "counts"), "CsparseMatrix")
  fh <- file.path(exp_dir, paste0(run, ".h5ad")); t0 <- Sys.time()
  write_h5ad(fh, cm, s@meta.data[, keep, drop = FALSE],
             if (length(rd)) SeuratObject::Embeddings(s, rd) else NULL)
  msg("%s  %.1f MB  (%.0f s)", basename(fh), file.size(fh)/1e6,
      as.numeric(Sys.time() - t0, units = "secs"))
}
