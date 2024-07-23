#' loadDatabase
#'
#' @return JASPAR预测所需数据库
#' @export
loadDatabase <- function() {
  suppressMessages(library(JASPAR2024))
  suppressMessages(library(TFBSTools))
  
  JASPAR2024 <- JASPAR2024()
  opts <- list()
  opts[["species"]] <- "Homo sapiens"
  opts[["matrixtype"]] <- "PWM"
  pwms <<- getMatrixSet(
    db(JASPAR2024),
    opts
  )
  
  opts[["matrixtype"]] <- "ICM"
  icms <<- getMatrixSet(
    db(JASPAR2024),
    opts
  )
}