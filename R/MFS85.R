#' @name MFS85
#' @title ALL Metabolic Subtypes
#' @description This R software package used GSVA to calculate the enrichment abundance of 85 metabolic pathways through transcriptome data of new samples, and based on the metabolic-feature-based subtypes (MFSs) and enrichment abundance of 85 metabolic pathways of the built-in training set, the new samples were classified as MFS1, or MFS2 or MFS3.
#' @param expr.data data.frame. Genes expression data.
#' @param gsva.methods character. The method of GSVA to estimate gene-set enrichment scores per sample. See \link{gsva}
#' @param ncores numeric. Number of threads used for GSVA analysis. See \link{gsva}
#' @param scale logical. Whether to scale the genes expression data.
#'
#' @import GSVA
#' @import pamr
#' @return result
#'
#'
#' @export
NULL
MFS85 <- function(expr.data,
                  gsva.methods = "gsva",
                  ncores = 1,
                  scale = FALSE
){
  #GSVA分析
  message("Step1: Now we will perform GSVA analysis... \n Please wait a minite!")
  #data(genelist)
  gsvaresult<- gsva(expr = as.matrix(expr.data),
                     gset.idx.list = genelist,
                     method=gsva.methods,
                     verbose=T,
                     parallel.sz=ncores
  )

  #scale数据
  if(scale){
    message("Step2: Scale expression data for pamr...")
    names <- dimnames(expr.data)
    expr.data <- t(scale(t(expr.data)))
    dimnames(expr.data) <- names
  }else{
    message("Step2: Don't scale expression data...")
  }

  #pamr分型
  message("Step3: classify the ALL metabolism type...")
  set.seed(10086)
  #data(train.exp)
  #data(train.class)
  x <- list(train.Exp=train.exp,train.MPS=train.class)
  data.pamr <- list(x=x$train.Exp,
                    y=x$train.MPS$cluster,
                    genenames=rownames(x$train.Exp)
  )
  model.train <- pamr.train(data.pamr)
  cv <- pamr.cv(model.train, data=data.pamr)
  thr <- cv$threshold[which.min(cv$error)]
  new.data <- gsvaresult[rownames(train.exp),]
  Pred <- pamr.predict(model.train, newx=new.data,
                       threshold=thr, type="class")
  Pred <- as.data.frame(Pred)
  rownames(Pred) <- colnames(new.data)
  result <- list(predict.class = Pred, gsvaresult = gsvaresult)
  return(result)
  message("MFS85 finished")
}

