#' Title
#'
#' @param rsID SNP的rs号
#' @param chr SNP所在染色体号
#' @param position SNP所在的碱基位置
#' @param ALT 风险等位基因
#' @param REF 参考等位基因
#' @param window SNP上下游扩展序列的窗口大小
#' @param strand 基因所在链的正反
#' @param threshold JASPAR预测显著结合的阈值
#' @param plotLogo 是否绘制转录因子Motif的LOGO图
#'
#' @return
#' @export
#'
#' @examples
#' singleJasparPredict("rs12345", 1, 1234567, "A", "T", 500, "+", 0.5, TRUE)
#' @export
singleJasparPredict <- function(rsID, chr, position, ALT, REF, window, strand, threshold, plotLogo) {

  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  suppressMessages(library(JASPAR2024))
  suppressMessages(library(TFBSTools))
  suppressMessages(library(dplyr))

  JASPAR2024 <- JASPAR2024()
  opts <- list()
  opts[["species"]] <- "Homo sapiens"
  opts[["matrixtype"]] <- "PWM"
  pwms <- getMatrixSet(
    db(JASPAR2024),
    opts
  )

  opts[["matrixtype"]] <- "ICM"
  icms <- getMatrixSet(
    db(JASPAR2024),
    opts
  )

  UPseq <- toString(getSeq(Hsapiens, paste0("chr", chr), start = position - window, end = position - 1))
  DOWNseq <- toString(getSeq(Hsapiens, paste0("chr", chr), start = position + 1, end = position + window))
  ALTseq <- paste0(UPseq, ALT, DOWNseq)
  REFseq <- paste0(UPseq, REF, DOWNseq)

  ALTres <- searchSeq(
    pwms, ALTseq, seqname=paste0(rsID, "[", ALT, "]"),
    min.score="0%", strand="+") %>%
    as.data.frame() %>%
    dplyr::select(ID, TF, class, start, end, ALTseqnames = seqnames, ALTabsScore = absScore, ALTrelScore = relScore, ALTsiteSeqs = siteSeqs)

  REFres <- searchSeq(
    pwms, REFseq, seqname=paste0(rsID, "[", REF, "]"),
    min.score="0%", strand="+") %>%
    as.data.frame() %>%
    dplyr::select(ID, TF, class, start, end, REFseqnames = seqnames, REFabsScore = absScore, REFrelScore = relScore, REFsiteSeqs = siteSeqs)

  outPut <- ALTres %>%
    dplyr::inner_join(REFres, by = c("ID", "TF", "class", "start", "end")) %>%
    mutate(class = case_when(
      ALTrelScore>=threshold & REFrelScore<threshold ~ "ALT Strong",
      ALTrelScore<threshold & REFrelScore>=threshold ~ "REF Strong",
      .default = "No Difference"
    )) %>%
    filter(class != "No Difference")

  if (plotLogo) {
    TFlist <- outPut %>% dplyr::distinct(ID, .keep_all = T)
    for (i in 1:nrow(TFlist)) {
      pdf(paste0(gsub("::", "-", TFlist$TF[i]), ".",  TFlist$ID[i], ".seqLogo.pdf"))
      seqLogo(icms[[which(names(icms) == TFlist$ID[i])]])
      dev.off()
    }
  }

  return(outPut)

}
