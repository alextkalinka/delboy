#' process_decision_boundary
#' 
#' Extract and smooth a decision boundary along the axes of fold-change and expression such that fold change is a convex, monotonically decreasing function of expression level.
#' 
#' @param data.grid A grid of predicted true and false positives along the plane formed by absolute log fold change and log10 expression.
#' 
#' @return A data frame.
#' @export
#' @importFrom grDevices chull
#' @importFrom smoothr smooth_ksmooth
process_decision_boundary <- function(data.grid){
  tryCatch({
    # Extract DB.
    db <- delboy::extract_grid_decision_boundary(data.grid)
    # Convex points only.
    db.chull <- grDevices::chull(as.matrix(db))
    # Densify and smooth.
    if(length(db.chull) > 5){
      db.ks <- smoothr::smooth_ksmooth(as.matrix(db[db.chull,]))
      db.chull <- grDevices::chull(db.ks)
      db_sm <- as.data.frame(db.ks[db.chull,])
    }else{
      db_sm <- db
    }
    colnames(db_sm) <- colnames(db)
    # Ensure DB is monotonically decreasing function of expression.
    db_sm <- delboy::smooth_decision_boundary(db_sm)
    db.chull <- grDevices::chull(as.matrix(db_sm))
    db_sm <- db_sm[db.chull,]
    colnames(db_sm) <- colnames(db)
    return(db_sm)
  },
  error = function(e) stop(paste("unable to process decision boundary:",e))
  )
}
