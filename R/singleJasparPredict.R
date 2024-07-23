#' single Jaspar Predict
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
#' @return 等位基因差异结合的转录因子列表
#' @export
singleJasparPredict <- function(rsID, chr, position, ALT, REF, window, strand, threshold, plotLogo = TRUE) {

  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  suppressMessages(library(TFBSTools))
  suppressMessages(library(dplyr))
  
  if ((!"icms"%in%ls(name = ".GlobalEnv")) | (!"pwms"%in%ls(name = ".GlobalEnv"))) {
    print("Please run loadDatabase() first!")
  }
  
  SNPseq <- toString(getSeq(Hsapiens, paste0("chr", chr), start = position, end = position))
  if (SNPseq != REF) {
    print("Please check reference allele and alternative allele!")
  }
    
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
