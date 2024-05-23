#' Best combination of normalization and imputation method
#'
#' @description This function will provide the best combinations of normalization
#' and imputation methods for the user given dataset based on the intragroup
#' variation evaluation parameters called PCV, PEV and PMAD.
#'
#' @param data_input Label-free proteomics expression data as a dataframe
#' @param groups Group information about the input data
#'
#' @details Label-free LC-MS proteomics expression data is often affected by heterogeneity and missing values. 
#' Normalization and missing value imputation are the commonly used techniques to solve these issues and make the dataset suitable for further downstream analysis. 
#' This function provides the best combination of normalization and imputation methods for the dataset, choosing from the three normalization methods (vsn, loess, and rlr) and three imputation methods (knn, lls, svd). 
#' The intragroup variation evaluation measures named pooled co-efficient of variance (PCV), pooled estimate of variance (PEV) and pooled median absolute deviation (PMAD) are used for selecting the best combination of normalization and imputation method for the given dataset.
#' It will return the best combinations based on each evaluation parameters of
#' PCV, PEV, and PMAD.
#'
#' Along with this, the user can get all three normalized datasets, nine combinations of normalized and missing values imputed datasets, and the PCV, PEV, and PMAD result values.
#'
#' @returns
#'  This function gives the list  which consist of following results.
#'
#' `Best Combinations`  The best combinations based on each PCV, PEV and PMAD
#' for the given dataset.
#'
#' `PCV Result` Values of groupwise PCV, overall PCV, PCV mean, PCV median and
#'              PCV standard deviation for all combinations.
#'
#' `PEV Result` Values of groupwise PEV, overall PEV, PEV mean, PEV median and
#'              PEV standard deviation for all combinations.
#'
#' `PMAD Result` Values of groupwise PMAD, overall PMAD, PMAD mean, PMAD median
#'               and PMAD standard deviation for all combinations.
#'
#' `vsn_data` The `vsn` normalized dataset
#'
#' `loess_data` The `loess` normalized dataset
#'
#' `rlr_data` The `rlr` normalized dataset
#'
#' `knn_vsn_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `knn` method.
#'
#' `knn_loess_data` The dataset normalized by `loess` method and missing values imputed
#'              by `knn` method.
#'
#' `knn_rlr_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `knn` method.
#'
#' `lls_vsn_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `lls` method.
#'
#' `lls_loess_data` The dataset normalized by `loess` method and missing values imputed
#'              by `lls` method.
#'
#' `lls_rlr_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `lls` method.
#'
#' `svd_vsn_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `svd` method.
#'
#' `svd_loess_data` The dataset normalized by `loess` method and missing values imputed
#'              by `svd` method.
#'
#' `svd_rlr_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `svd` method.
#'
#' @author
#'  Dr Sudhir Srivastava ("Sudhir.Srivastava@icar.gov.in")
#'
#'  Kabilan S ("kabilan151414@gmail.com")
#'
#' @export
#'
#' @examples
#' \donttest{
#' result <- best_combination(yeast_data, yeast_groups)
#' result$`Best combinations`
#' result$`PCV Result`
#' result$`PMAD Result`
#' result$`knn_rlr_data`
#' }

#Main function for finding out the top three combinations
best_combination <- function (data_input, groups){
  
  sink("output.txt")

  Type <- Group <- name <- value <- . <- value_Mean <- median <-
    value_Median <- sd <- value_SD <- PCV_mean <- PCV_median <- PCV_sd <-
    PEV_mean <- PEV_median <- PEV_sd <- PMAD_mean <- PMAD_median <- PMAD_sd <-NULL
  
  #Converting all zeros to NAs
  data_input[data_input == 0] <- NA
  
  #Complete data function (To remove complete missing values in the group)
  complete_data_fn <- function (data, groups){
    
    rowid <- id <- group <- mass <- n <- name <- NULL
    
    #Rename the groups data
    new_colnames <- c("name", "group")
    colnames(groups)<-new_colnames
    
    #Rename the first column of data_input
    data_input1 <- data
    colnames(data_input1)[1] <- "rowid"
    
    #Removing groupwise missing rows in a dataframe
    com_data <- data_input1 %>%
      dplyr::group_by(rowid) %>%
      dplyr::mutate(id = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(-c(rowid, id), values_to = "mass") %>%
      dplyr::inner_join(groups, by = "name") %>%
      dplyr::add_count(rowid, id, group, wt = !is.na(mass)) %>%
      dplyr::group_by(rowid, id) %>%
      dplyr::filter(!any(n == 0)) %>%
      dplyr::ungroup() %>%
      dplyr::select(!c(group, n)) %>%
      tidyr::pivot_wider(names_from = name, values_from = mass) %>%
      dplyr::select(-id)
    
    return(com_data)
  }
  
  #Getting the complete data
  com_data <- complete_data_fn (data_input, groups)
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  
  #Grouping of dataframe
  grouping_data<-function(df){                      #df= dataframe
    df_col<-ncol(df)                                #calculates no. of columns in dataframe
    x <- table(groups$Groups)               #Extract the unique group value
    s <- unique(x)
    groups<-sort(rep(0:((df_col/s)-1),s))
    id<-list()                                      #creates empty list
    for (i in 1:length(unique(groups))){
      id[[i]]<-which(groups == unique(groups)[i])}  #creates list of groups
    names(id)<-paste0("id",unique(groups))          #assigns group based names to the list "id"
    data<-list()                                    #creates empty list
    for (i in 1:length(id)){
      data[[i]]<-df[,id[[i]]]}                      #creates list of dataframe columns sorted by groups
    names(data)<-paste0("data",unique(groups))      #assigns group based names to the list "data"
    return(data)}
  
  #Grouping of dataframe as a triplicate groups
  group_data <- grouping_data(com_data2)
  
  #VSN Normalization function
  VSN_Norm <- function(dat) {
    dat<-as.data.frame(dat)
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
  }
  
  #VSN normalized data
  vsn.dat <- do.call("cbind", lapply(group_data, VSN_Norm))
  
  #Loess normalization
  LOESS_Norm <- function(dat) {
    newdata<- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method="fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  
  #Loess normalized data
  loess.dat <- do.call("cbind", lapply(group_data, LOESS_Norm))
  
  #RLR normalization
  RLR_Norm <- function (dat) {
    log2Matrix <- log2(dat)
    log2Matrix <- as.matrix(log2Matrix)
    
    # Extract the numeric component from the list
    log2Matrix <- unlist(log2Matrix)
    log2Matrix[is.infinite(log2Matrix)] <- 0
    
    sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
    calculateRLMForCol <- function(colIndex, sampleLog2Median,
                                   log2Matrix) {
      lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~
                           sampleLog2Median, na.action = stats::na.exclude, maxit = 20)
      coeffs <- lrFit$coefficients
      coefIntercept <- coeffs[1]
      coefSlope <- coeffs[2]
      globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept)/coefSlope
      globalFittedRLRCol
    }
    
    data <- log2Matrix
    globalFittedRLR <- vapply(seq_len(ncol(data)), calculateRLMForCol,
                              rep(0, nrow(data)), sampleLog2Median = sampleLog2Median,
                              log2Matrix = data)
    colnames(globalFittedRLR) <- colnames(dat)
    globalFittedRLR[globalFittedRLR < 0] <- 0
    globalFittedRLR
  }
  
  #RLR normalized data
  rlr.dat <- do.call("cbind", lapply(group_data, RLR_Norm))
  
  #Grouping of normalized datasets
  vsn_group_data <- grouping_data(vsn.dat)
  
  loess_group_data <- grouping_data(loess.dat)
  
  rlr_group_data <- grouping_data(rlr.dat)
  
  #Imputation of normalized datasets
  #KNN imputation
  KNN_Imputation <- function (dat)
  {
    resultkNN <- VIM::kNN(dat, numFun = laeken::weightedMean, weightDist = TRUE,
                          imp_var = FALSE, k= 10)
    return(resultkNN)
  }
  
  #LLS imputation
  LLS_Imputation <- function (dat)
  {
    resultLLS <- pcaMethods::llsImpute(dat, k=2, correlation = "pearson", allVariables = TRUE)
    dataSet.imputed <- resultLLS@completeObs
    return(dataSet.imputed)
  }
  
  #SVD imputation
  SVD_Imputation <- function (dat)
  {
    resultSVD <- pcaMethods::pca(dat, method = "svdImpute", nPcs = 2)
    dataSet.imputed <- resultSVD@completeObs
    return(dataSet.imputed)
  }
  
  #KNN and VSN
  knn.vsn.dat <- do.call("cbind", lapply(vsn_group_data, KNN_Imputation))
  
  #LLS and VSN
  lls.vsn.dat <- do.call("cbind", lapply(vsn_group_data, LLS_Imputation))
  
  #SVD and VSN
  svd.vsn.dat <- do.call("cbind", lapply(vsn_group_data, SVD_Imputation))
  
  #KNN and LOESS
  knn.loess.dat <- do.call("cbind", lapply(loess_group_data, KNN_Imputation))
  
  #LLS and LOESS
  lls.loess.dat <- do.call("cbind", lapply(loess_group_data, LLS_Imputation))
  
  #SVD and LOESS
  svd.loess.dat <- do.call("cbind", lapply(loess_group_data, SVD_Imputation))
  
  #KNN and RLR
  knn.rlr.dat <- do.call("cbind", lapply(rlr_group_data, KNN_Imputation))
  
  #LLS and RLR
  lls.rlr.dat <- do.call("cbind", lapply(rlr_group_data, LLS_Imputation))
  
  #SVD and RLR
  svd.rlr.dat <- do.call("cbind", lapply(rlr_group_data, SVD_Imputation))
  
  #Transposing the sample to wide format
  new_sample <- as.data.frame(t(groups))
  names(new_sample) <- new_sample[1,]
  sample <- new_sample[-1,]
  
  #Changing all data files column names
  new_colnames <- colnames(sample)
  colnames(knn.vsn.dat) <- new_colnames
  colnames(knn.loess.dat) <- new_colnames
  colnames(knn.rlr.dat) <- new_colnames
  colnames(lls.vsn.dat) <- new_colnames
  colnames(lls.loess.dat) <- new_colnames
  colnames(lls.rlr.dat) <- new_colnames
  colnames(svd.vsn.dat) <- new_colnames
  colnames(svd.loess.dat) <- new_colnames
  colnames(svd.rlr.dat) <- new_colnames
  
  #Group wise PCV
  Group_data_PCV = function(data, groups){
    PCV = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      CVs = apply(data, 1, sd, na.rm = FALSE)/
        rowMeans(tempData, na.rm = FALSE)
      PCV_mean[group] = mean(CVs, na.rm = T)
      PCV_median[group] = median(CVs, na.rm = T)
      PCV_sd[group] = sd(CVs, na.rm = T)
    }
    Group_data_PCV_list = list("Group_data_PCV_mean" = as.data.frame(PCV_mean),
                               "Group_data_PCV_median" = as.data.frame(PCV_median),
                               "Group_data_PCV_sd" = as.data.frame(PCV_sd))
    return(Group_data_PCV_list)
  }
  
  #PCV calculation_groupwise
  test1 <- Group_data_PCV(knn.vsn.dat, sample)
  knn_vsn_PCV_mean <- cbind(Type ="knn_vsn", test1$Group_data_PCV_mean)
  knn_vsn_PCV_median <- cbind(Type ="knn_vsn", test1$Group_data_PCV_median)
  knn_vsn_PCV_sd <- cbind(Type ="knn_vsn", test1$Group_data_PCV_sd)
  
  test2 <- Group_data_PCV(knn.loess.dat, sample)
  knn_loess_PCV_mean <- cbind(Type ="knn_loess", test2$Group_data_PCV_mean)
  knn_loess_PCV_median <- cbind(Type ="knn_loess", test2$Group_data_PCV_median)
  knn_loess_PCV_sd <- cbind(Type ="knn_loess", test2$Group_data_PCV_sd)
  
  test3 <- Group_data_PCV(knn.rlr.dat, sample)
  knn_rlr_PCV_mean <- cbind(Type ="knn_rlr", test3$Group_data_PCV_mean)
  knn_rlr_PCV_median <- cbind(Type ="knn_rlr", test3$Group_data_PCV_median)
  knn_rlr_PCV_sd <- cbind(Type ="knn_rlr", test3$Group_data_PCV_sd)
  
  test4 <- Group_data_PCV(lls.vsn.dat, sample)
  lls_vsn_PCV_mean <- cbind(Type ="lls_vsn", test4$Group_data_PCV_mean)
  lls_vsn_PCV_median <- cbind(Type ="lls_vsn", test4$Group_data_PCV_median)
  lls_vsn_PCV_sd <- cbind(Type ="lls_vsn", test4$Group_data_PCV_sd)
  
  test5 <- Group_data_PCV(lls.loess.dat, sample)
  lls_loess_PCV_mean <- cbind(Type ="lls_loess", test5$Group_data_PCV_mean)
  lls_loess_PCV_median <- cbind(Type ="lls_loess", test5$Group_data_PCV_median)
  lls_loess_PCV_sd <- cbind(Type ="lls_loess", test5$Group_data_PCV_sd)
  
  test6 <- Group_data_PCV(lls.rlr.dat, sample)
  lls_rlr_PCV_mean <- cbind(Type ="lls_rlr", test6$Group_data_PCV_mean)
  lls_rlr_PCV_median <- cbind(Type ="lls_rlr", test6$Group_data_PCV_median)
  lls_rlr_PCV_sd <- cbind(Type ="lls_rlr", test6$Group_data_PCV_sd)
  
  test7 <- Group_data_PCV(svd.vsn.dat, sample)
  svd_vsn_PCV_mean <- cbind(Type ="svd_vsn", test7$Group_data_PCV_mean)
  svd_vsn_PCV_median <- cbind(Type ="svd_vsn", test7$Group_data_PCV_median)
  svd_vsn_PCV_sd <- cbind(Type ="svd_vsn", test7$Group_data_PCV_sd)
  
  test8 <- Group_data_PCV(svd.loess.dat, sample)
  svd_loess_PCV_mean <- cbind(Type ="svd_loess", test8$Group_data_PCV_mean)
  svd_loess_PCV_median <- cbind(Type ="svd_loess", test8$Group_data_PCV_median)
  svd_loess_PCV_sd <- cbind(Type ="svd_loess", test8$Group_data_PCV_sd)
  
  test9 <- Group_data_PCV(svd.rlr.dat, sample)
  svd_rlr_PCV_mean <- cbind(Type ="svd_rlr", test9$Group_data_PCV_mean)
  svd_rlr_PCV_median <- cbind(Type ="svd_rlr", test9$Group_data_PCV_median)
  svd_rlr_PCV_sd <- cbind(Type ="svd_rlr", test9$Group_data_PCV_sd)
  
  ###PCV_mean
  #Combining all the above results
  total_pcv_mean <- plyr::rbind.fill(knn_vsn_PCV_mean, knn_loess_PCV_mean, knn_rlr_PCV_mean, 
                                     lls_vsn_PCV_mean, lls_loess_PCV_mean, lls_rlr_PCV_mean, 
                                     svd_vsn_PCV_mean, svd_loess_PCV_mean, svd_rlr_PCV_mean)
  
  #Separating the results groupwise
  total_Group_data_PCV_mean <- total_pcv_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_mean)
  
  #Extract the top combination in each group
  total_Group_data_PCV_mean2<-
    total_Group_data_PCV_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PCV_mean <- subset(total_Group_data_PCV_mean2, select = -row)
  
  ###PCV_median
  #Combining all the above results
  total_pcv_median <- plyr::rbind.fill(knn_vsn_PCV_median, knn_loess_PCV_median, knn_rlr_PCV_median, 
                                       lls_vsn_PCV_median, lls_loess_PCV_median, lls_rlr_PCV_median, 
                                       svd_vsn_PCV_median, svd_loess_PCV_median, svd_rlr_PCV_median)
  
  #Separating the results groupwise
  total_Group_data_PCV_median <- total_pcv_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_median)
  
  #Extract the top combination in each group
  total_Group_data_PCV_median2<-
    total_Group_data_PCV_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result 
  final_Group_data_PCV_median <- subset(total_Group_data_PCV_median2, select = -row)
  
  ###PCV_sd
  #Combining all the above results
  total_pcv_sd <- plyr::rbind.fill(knn_vsn_PCV_sd, knn_loess_PCV_sd, knn_rlr_PCV_sd, 
                                   lls_vsn_PCV_sd, lls_loess_PCV_sd, lls_rlr_PCV_sd, 
                                   svd_vsn_PCV_sd, svd_loess_PCV_sd, svd_rlr_PCV_sd)
  
  #Separating the results groupwise
  total_Group_data_PCV_sd <- total_pcv_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_sd)
  
  #Extract the top combination in each group
  total_Group_data_PCV_sd2<-
    total_Group_data_PCV_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PCV_sd <- subset(total_Group_data_PCV_sd2, select = -row)
  
  #Overall data PCV
  Total_data_PCV = function(data){
    CVs = apply(data, 1, sd, na.rm = FALSE)/
      rowMeans(data, na.rm = FALSE)
    PCV_mean = mean(CVs, na.rm = T)
    PCV_median =  median(CVs, na.rm = T)
    PCV_sd = sd(CVs, na.rm = T)
    PCV_list = list("PCV_mean" = as.data.frame(PCV_mean),
                    "PCV_median" = as.data.frame(PCV_median),
                    "PCV_sd" = as.data.frame(PCV_sd))
    return(PCV_list)
  }
  
  
  #Overall PCV_mean estimation
  test1 <- Total_data_PCV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PCV_mean)
  
  test2 <- Total_data_PCV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PCV_mean)
  
  test3 <- Total_data_PCV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PCV_mean)
  
  test4 <- Total_data_PCV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PCV_mean)
  
  test5 <- Total_data_PCV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PCV_mean)
  
  test6 <- Total_data_PCV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PCV_mean)
  
  test7 <- Total_data_PCV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PCV_mean)
  
  test8 <- Total_data_PCV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PCV_mean)
  
  test9 <- Total_data_PCV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PCV_mean)
  
  #Combining all the above results
  total_pcv_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_mean <-total_pcv_overall_mean2%>%dplyr::slice_min(PCV_mean, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PCV_mean", "Overall_value.PCV_mean")
  colnames(total_pcv_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PCV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PCV_median)
  
  test2 <- Total_data_PCV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PCV_median)
  
  test3 <- Total_data_PCV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PCV_median)
  
  test4 <- Total_data_PCV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PCV_median)
  
  test5 <- Total_data_PCV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PCV_median)
  
  test6 <- Total_data_PCV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PCV_median)
  
  test7 <- Total_data_PCV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PCV_median)
  
  test8 <- Total_data_PCV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PCV_median)
  
  test9 <- Total_data_PCV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PCV_median)
  
  #Combining all the above results
  total_pcv_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_median <-total_pcv_overall_median2%>%dplyr::slice_min(PCV_median, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PCV_median", "Overall_value.PCV_median")
  colnames(total_pcv_overall_median) <- original_cols
  
  #Overall PCV_sd estimation
  test1 <- Total_data_PCV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PCV_sd)
  
  test2 <- Total_data_PCV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PCV_sd)
  
  test3 <- Total_data_PCV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PCV_sd)
  
  test4 <- Total_data_PCV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PCV_sd)
  
  test5 <- Total_data_PCV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PCV_sd)
  
  test6 <- Total_data_PCV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PCV_sd)
  
  test7 <- Total_data_PCV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PCV_sd)
  
  test8 <- Total_data_PCV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PCV_sd)
  
  test9 <- Total_data_PCV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PCV_sd)
  
  #Combining all the above results
  total_pcv_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_sd <-total_pcv_overall_sd2%>%dplyr::slice_min(PCV_sd, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PCV_sd", "Overall_value.PCV_sd")
  colnames(total_pcv_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PCV<- cbind(final_Group_data_PCV_mean, final_Group_data_PCV_median,
                     final_Group_data_PCV_sd, total_pcv_overall_mean,
                     total_pcv_overall_median, total_pcv_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PCV), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PCV_names1 <- result_PCV[, -col_indices]
  n <- matrix(t(result_PCV_names1), ncol=1)
  all_elements <- unlist(n)
  
  # Get the frequency of each element in the dataframe
  freq <- table(all_elements)
  
  # Find the most occurring element
  PCV_best_combination <- names(freq)[which.max(freq)]
  
  #Groupwise PEV estimation
  Group_data_PEV = function(data1, groups){
    data <- as.matrix(data1)
    PEV = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      rowNonNACnt = rowSums(!is.na(tempData)) - 1
      EV = rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
      PEV_mean[group] = mean(EV, na.rm = T)
      PEV_median[group] = median(EV, na.rm = T)
      PEV_sd[group] = sd(EV, na.rm = T)
    }
    Group_data_PEV_list = list("Group_data_PEV_mean" = as.data.frame(PEV_mean),
                               "Group_data_PEV_median" = as.data.frame(PEV_median),
                               "Group_data_PEV_sd" = as.data.frame(PEV_sd))
    return(Group_data_PEV_list)
  }
  
  
  #PEV calculation_groupwise
  test1 <- Group_data_PEV(knn.vsn.dat, sample)
  knn_vsn_PEV_mean <- cbind(Type ="knn_vsn", test1$Group_data_PEV_mean)
  knn_vsn_PEV_median <- cbind(Type ="knn_vsn", test1$Group_data_PEV_median)
  knn_vsn_PEV_sd <- cbind(Type ="knn_vsn", test1$Group_data_PEV_sd)
  
  test2 <- Group_data_PEV(knn.loess.dat, sample)
  knn_loess_PEV_mean <- cbind(Type ="knn_loess", test2$Group_data_PEV_mean)
  knn_loess_PEV_median <- cbind(Type ="knn_loess", test2$Group_data_PEV_median)
  knn_loess_PEV_sd <- cbind(Type ="knn_loess", test2$Group_data_PEV_sd)
  
  test3 <- Group_data_PEV(knn.rlr.dat, sample)
  knn_rlr_PEV_mean <- cbind(Type ="knn_rlr", test3$Group_data_PEV_mean)
  knn_rlr_PEV_median <- cbind(Type ="knn_rlr", test3$Group_data_PEV_median)
  knn_rlr_PEV_sd <- cbind(Type ="knn_rlr", test3$Group_data_PEV_sd)
  
  test4 <- Group_data_PEV(lls.vsn.dat, sample)
  lls_vsn_PEV_mean <- cbind(Type ="lls_vsn", test4$Group_data_PEV_mean)
  lls_vsn_PEV_median <- cbind(Type ="lls_vsn", test4$Group_data_PEV_median)
  lls_vsn_PEV_sd <- cbind(Type ="lls_vsn", test4$Group_data_PEV_sd)
  
  test5 <- Group_data_PEV(lls.loess.dat, sample)
  lls_loess_PEV_mean <- cbind(Type ="lls_loess", test5$Group_data_PEV_mean)
  lls_loess_PEV_median <- cbind(Type ="lls_loess", test5$Group_data_PEV_median)
  lls_loess_PEV_sd <- cbind(Type ="lls_loess", test5$Group_data_PEV_sd)
  
  test6 <- Group_data_PEV(lls.rlr.dat, sample)
  lls_rlr_PEV_mean <- cbind(Type ="lls_rlr", test6$Group_data_PEV_mean)
  lls_rlr_PEV_median <- cbind(Type ="lls_rlr", test6$Group_data_PEV_median)
  lls_rlr_PEV_sd <- cbind(Type ="lls_rlr", test6$Group_data_PEV_sd)
  
  test7 <- Group_data_PEV(svd.vsn.dat, sample)
  svd_vsn_PEV_mean <- cbind(Type ="svd_vsn", test7$Group_data_PEV_mean)
  svd_vsn_PEV_median <- cbind(Type ="svd_vsn", test7$Group_data_PEV_median)
  svd_vsn_PEV_sd <- cbind(Type ="svd_vsn", test7$Group_data_PEV_sd)
  
  test8 <- Group_data_PEV(svd.loess.dat, sample)
  svd_loess_PEV_mean <- cbind(Type ="svd_loess", test8$Group_data_PEV_mean)
  svd_loess_PEV_median <- cbind(Type ="svd_loess", test8$Group_data_PEV_median)
  svd_loess_PEV_sd <- cbind(Type ="svd_loess", test8$Group_data_PEV_sd)
  
  test9 <- Group_data_PEV(svd.rlr.dat, sample)
  svd_rlr_PEV_mean <- cbind(Type ="svd_rlr", test9$Group_data_PEV_mean)
  svd_rlr_PEV_median <- cbind(Type ="svd_rlr", test9$Group_data_PEV_median)
  svd_rlr_PEV_sd <- cbind(Type ="svd_rlr", test9$Group_data_PEV_sd)
  
  ###PEV_mean
  #Combining all the above results
  total_pev_mean <- plyr::rbind.fill(knn_vsn_PEV_mean, knn_loess_PEV_mean, knn_rlr_PEV_mean, 
                                     lls_vsn_PEV_mean, lls_loess_PEV_mean, lls_rlr_PEV_mean, 
                                     svd_vsn_PEV_mean, svd_loess_PEV_mean, svd_rlr_PEV_mean)
  
  #Separating the results groupwise
  total_Group_data_PEV_mean <- total_pev_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_mean)
  
  #Extract the top combination in each group
  total_Group_data_PEV_mean2<-
    total_Group_data_PEV_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PEV_mean <- subset(total_Group_data_PEV_mean2, select = -row)
  
  ###PEV_median
  #Combining all the above results
  total_pev_median <- plyr::rbind.fill(knn_vsn_PEV_median, knn_loess_PEV_median, knn_rlr_PEV_median, 
                                       lls_vsn_PEV_median, lls_loess_PEV_median, lls_rlr_PEV_median, 
                                       svd_vsn_PEV_median, svd_loess_PEV_median, svd_rlr_PEV_median)
  
  #Separating the results groupwise
  total_Group_data_PEV_median <- total_pev_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_median)
  
  #Extract the top combination in each group
  total_Group_data_PEV_median2<-
    total_Group_data_PEV_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result 
  final_Group_data_PEV_median <- subset(total_Group_data_PEV_median2, select = -row)
  
  ###PEV_sd
  #Combining all the above results
  total_pev_sd <- plyr::rbind.fill(knn_vsn_PEV_sd, knn_loess_PEV_sd, knn_rlr_PEV_sd, 
                                   lls_vsn_PEV_sd, lls_loess_PEV_sd, lls_rlr_PEV_sd, 
                                   svd_vsn_PEV_sd, svd_loess_PEV_sd, svd_rlr_PEV_sd)
  
  #Separating the results groupwise
  total_Group_data_PEV_sd <- total_pev_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_sd)
  
  #Extract the top combination in each group
  total_Group_data_PEV_sd2<-
    total_Group_data_PEV_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PEV_sd <- subset(total_Group_data_PEV_sd2, select = -row)
  
  #Overall PEV function
  Total_data_PEV = function(data1){
    data <- as.matrix(data1)
    tempData = data[]
    rowNonNACnt = rowSums(!is.na(tempData)) - 1
    EVs = rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
    PEV_mean = mean(EVs, na.rm = FALSE)
    PEV_median = median(EVs, na.rm = FALSE)
    PEV_sd = sd(EVs, na.rm = FALSE)
    PEV_list = list("PEV_mean" = as.data.frame(PEV_mean),
                    "PEV_median" = as.data.frame(PEV_median),
                    "PEV_sd" = as.data.frame(PEV_sd))
    return(PEV_list)
  }
  
  #Overall PEV_mean estimation
  test1 <- Total_data_PEV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PEV_mean)
  
  test2 <- Total_data_PEV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PEV_mean)
  
  test3 <- Total_data_PEV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PEV_mean)
  
  test4 <- Total_data_PEV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PEV_mean)
  
  test5 <- Total_data_PEV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PEV_mean)
  
  test6 <- Total_data_PEV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PEV_mean)
  
  test7 <- Total_data_PEV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PEV_mean)
  
  test8 <- Total_data_PEV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PEV_mean)
  
  test9 <- Total_data_PEV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PEV_mean)
  
  
  #Combining all the above results
  total_pev_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_mean <-total_pev_overall_mean2%>%dplyr::slice_min(PEV_mean, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PEV_mean", "Overall_value.PEV_mean")
  colnames(total_pev_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PEV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PEV_median)
  
  test2 <- Total_data_PEV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PEV_median)
  
  test3 <- Total_data_PEV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PEV_median)
  
  test4 <- Total_data_PEV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PEV_median)
  
  test5 <- Total_data_PEV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PEV_median)
  
  test6 <- Total_data_PEV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PEV_median)
  
  test7 <- Total_data_PEV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PEV_median)
  
  test8 <- Total_data_PEV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PEV_median)
  
  test9 <- Total_data_PEV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PEV_median)
  
  #Combining all the above results
  total_pev_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_median <-total_pev_overall_median2%>%dplyr::slice_min(PEV_median, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PEV_median", "Overall_value.PEV_median")
  colnames(total_pev_overall_median) <- original_cols
  
  #Overall PEV_sd estimation
  test1 <- Total_data_PEV(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PEV_sd)
  
  test2 <- Total_data_PEV(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PEV_sd)
  
  test3 <- Total_data_PEV(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PEV_sd)
  
  test4 <- Total_data_PEV(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PEV_sd)
  
  test5 <- Total_data_PEV(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PEV_sd)
  
  test6 <- Total_data_PEV(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PEV_sd)
  
  test7 <- Total_data_PEV(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PEV_sd)
  
  test8 <- Total_data_PEV(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PEV_sd)
  
  test9 <- Total_data_PEV(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PEV_sd)
  
  #Combining all the above results
  total_pev_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_sd <-total_pev_overall_sd2%>%dplyr::slice_min(PEV_sd, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PEV_sd", "Overall_value.PEV_sd")
  colnames(total_pev_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PEV<- cbind(final_Group_data_PEV_mean, final_Group_data_PEV_median,
                     final_Group_data_PEV_sd, total_pev_overall_mean,
                     total_pev_overall_median, total_pev_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PEV), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PEV_names1 <- result_PEV[, -col_indices]
  n <- matrix(t(result_PEV_names1), ncol=1)
  all_elements <- unlist(n)
  
  # Get the frequency of each element in the dataframe
  freq <- table(all_elements)
  
  # Find the most occurring element
  PEV_best_combination <- names(freq)[which.max(freq)]
  
  #Groupwise PMAD function
  Group_data_PMAD = function(data1, groups){
    data <- as.matrix(data1)
    PMAD = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      MAD = matrixStats::rowMads(tempData, na.rm = FALSE)
      PMAD_mean[group] = mean(MAD, na.rm = T)
      PMAD_median[group] = median(MAD, na.rm = T)
      PMAD_sd[group] = sd(MAD, na.rm = T)
    }
    Group_data_PMAD_list = list("Group_data_PMAD_mean" = as.data.frame(PMAD_mean),
                                "Group_data_PMAD_median" = as.data.frame(PMAD_median),
                                "Group_data_PMAD_sd" = as.data.frame(PMAD_sd))
    return(Group_data_PMAD_list)
  }
  
  #PMAD calculation_groupwise
  test1 <- Group_data_PMAD(knn.vsn.dat, sample)
  knn_vsn_PMAD_mean <- cbind(Type ="knn_vsn", test1$Group_data_PMAD_mean)
  knn_vsn_PMAD_median <- cbind(Type ="knn_vsn", test1$Group_data_PMAD_median)
  knn_vsn_PMAD_sd <- cbind(Type ="knn_vsn", test1$Group_data_PMAD_sd)
  
  test2 <- Group_data_PMAD(knn.loess.dat, sample)
  knn_loess_PMAD_mean <- cbind(Type ="knn_loess", test2$Group_data_PMAD_mean)
  knn_loess_PMAD_median <- cbind(Type ="knn_loess", test2$Group_data_PMAD_median)
  knn_loess_PMAD_sd <- cbind(Type ="knn_loess", test2$Group_data_PMAD_sd)
  
  test3 <- Group_data_PMAD(knn.rlr.dat, sample)
  knn_rlr_PMAD_mean <- cbind(Type ="knn_rlr", test3$Group_data_PMAD_mean)
  knn_rlr_PMAD_median <- cbind(Type ="knn_rlr", test3$Group_data_PMAD_median)
  knn_rlr_PMAD_sd <- cbind(Type ="knn_rlr", test3$Group_data_PMAD_sd)
  
  test4 <- Group_data_PMAD(lls.vsn.dat, sample)
  lls_vsn_PMAD_mean <- cbind(Type ="lls_vsn", test4$Group_data_PMAD_mean)
  lls_vsn_PMAD_median <- cbind(Type ="lls_vsn", test4$Group_data_PMAD_median)
  lls_vsn_PMAD_sd <- cbind(Type ="lls_vsn", test4$Group_data_PMAD_sd)
  
  test5 <- Group_data_PMAD(lls.loess.dat, sample)
  lls_loess_PMAD_mean <- cbind(Type ="lls_loess", test5$Group_data_PMAD_mean)
  lls_loess_PMAD_median <- cbind(Type ="lls_loess", test5$Group_data_PMAD_median)
  lls_loess_PMAD_sd <- cbind(Type ="lls_loess", test5$Group_data_PMAD_sd)
  
  test6 <- Group_data_PMAD(lls.rlr.dat, sample)
  lls_rlr_PMAD_mean <- cbind(Type ="lls_rlr", test6$Group_data_PMAD_mean)
  lls_rlr_PMAD_median <- cbind(Type ="lls_rlr", test6$Group_data_PMAD_median)
  lls_rlr_PMAD_sd <- cbind(Type ="lls_rlr", test6$Group_data_PMAD_sd)
  
  test7 <- Group_data_PMAD(svd.vsn.dat, sample)
  svd_vsn_PMAD_mean <- cbind(Type ="svd_vsn", test7$Group_data_PMAD_mean)
  svd_vsn_PMAD_median <- cbind(Type ="svd_vsn", test7$Group_data_PMAD_median)
  svd_vsn_PMAD_sd <- cbind(Type ="svd_vsn", test7$Group_data_PMAD_sd)
  
  test8 <- Group_data_PMAD(svd.loess.dat, sample)
  svd_loess_PMAD_mean <- cbind(Type ="svd_loess", test8$Group_data_PMAD_mean)
  svd_loess_PMAD_median <- cbind(Type ="svd_loess", test8$Group_data_PMAD_median)
  svd_loess_PMAD_sd <- cbind(Type ="svd_loess", test8$Group_data_PMAD_sd)
  
  test9 <- Group_data_PMAD(svd.rlr.dat, sample)
  svd_rlr_PMAD_mean <- cbind(Type ="svd_rlr", test9$Group_data_PMAD_mean)
  svd_rlr_PMAD_median <- cbind(Type ="svd_rlr", test9$Group_data_PMAD_median)
  svd_rlr_PMAD_sd <- cbind(Type ="svd_rlr", test9$Group_data_PMAD_sd)
  
  ###PMAD_mean
  #Combining all the above results
  total_pmad_mean <- plyr::rbind.fill(knn_vsn_PMAD_mean, knn_loess_PMAD_mean, knn_rlr_PMAD_mean, 
                                      lls_vsn_PMAD_mean, lls_loess_PMAD_mean, lls_rlr_PMAD_mean, 
                                      svd_vsn_PMAD_mean, svd_loess_PMAD_mean, svd_rlr_PMAD_mean)
  
  #Separating the results groupwise
  total_Group_data_PMAD_mean <- total_pmad_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_mean)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_mean2<-
    total_Group_data_PMAD_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PMAD_mean <- subset(total_Group_data_PMAD_mean2, select = -row)
  
  ###PMAD_median
  #Combining all the above results
  total_pmad_median <- plyr::rbind.fill(knn_vsn_PMAD_median, knn_loess_PMAD_median, knn_rlr_PMAD_median, 
                                        lls_vsn_PMAD_median, lls_loess_PMAD_median, lls_rlr_PMAD_median, 
                                        svd_vsn_PMAD_median, svd_loess_PMAD_median, svd_rlr_PMAD_median)
  
  #Separating the results groupwise
  total_Group_data_PMAD_median <- total_pmad_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_median)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_median2<-
    total_Group_data_PMAD_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result 
  final_Group_data_PMAD_median <- subset(total_Group_data_PMAD_median2, select = -row)
  
  ###PMAD_sd
  #Combining all the above results
  total_pmad_sd <- plyr::rbind.fill(knn_vsn_PMAD_sd, knn_loess_PMAD_sd, knn_rlr_PMAD_sd, 
                                    lls_vsn_PMAD_sd, lls_loess_PMAD_sd, lls_rlr_PMAD_sd, 
                                    svd_vsn_PMAD_sd, svd_loess_PMAD_sd, svd_rlr_PMAD_sd)
  
  #Separating the results groupwise
  total_Group_data_PMAD_sd <- total_pmad_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_sd)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_sd2<-
    total_Group_data_PMAD_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PMAD_sd <- subset(total_Group_data_PMAD_sd2, select = -row)
  
  #Overall PMAD function
  Total_data_PMAD = function(data1){
    data <- as.matrix(data1)
    tempData = data[]
    MADs = matrixStats::rowMads(tempData, na.rm = FALSE)
    PMAD_mean = mean(MADs, na.rm = T)
    PMAD_median = median(MADs, na.rm = T)
    PMAD_sd = sd(MADs, na.rm = T)
    PMAD_list = list("PMAD_mean" = as.data.frame(PMAD_mean),
                     "PMAD_median" = as.data.frame(PMAD_median),
                     "PMAD_sd" = as.data.frame(PMAD_sd))
    return(PMAD_list)
  }
  
  #Overall PMAD_mean estimation
  test1 <- Total_data_PMAD(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PMAD_mean)
  
  test2 <- Total_data_PMAD(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PMAD_mean)
  
  test3 <- Total_data_PMAD(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PMAD_mean)
  
  test4 <- Total_data_PMAD(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PMAD_mean)
  
  test5 <- Total_data_PMAD(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PMAD_mean)
  
  test6 <- Total_data_PMAD(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PMAD_mean)
  
  test7 <- Total_data_PMAD(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PMAD_mean)
  
  test8 <- Total_data_PMAD(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PMAD_mean)
  
  test9 <- Total_data_PMAD(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PMAD_mean)
  
  
  #Combining all the above results
  total_pmad_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_mean <-total_pmad_overall_mean2%>%dplyr::slice_min(PMAD_mean, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PMAD_mean", "Overall_value.PMAD_mean")
  colnames(total_pmad_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PMAD(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PMAD_median)
  
  test2 <- Total_data_PMAD(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PMAD_median)
  
  test3 <- Total_data_PMAD(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PMAD_median)
  
  test4 <- Total_data_PMAD(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PMAD_median)
  
  test5 <- Total_data_PMAD(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PMAD_median)
  
  test6 <- Total_data_PMAD(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PMAD_median)
  
  test7 <- Total_data_PMAD(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PMAD_median)
  
  test8 <- Total_data_PMAD(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PMAD_median)
  
  test9 <- Total_data_PMAD(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PMAD_median)
  
  #Combining all the above results
  total_pmad_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_median <-total_pmad_overall_median2%>%dplyr::slice_min(PMAD_median, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PMAD_median", "Overall_value.PMAD_median")
  colnames(total_pmad_overall_median) <- original_cols
  
  #Overall PMAD_sd estimation
  test1 <- Total_data_PMAD(knn.vsn.dat)
  data1 <- cbind(Type ="knn_vsn", test1$PMAD_sd)
  
  test2 <- Total_data_PMAD(knn.loess.dat)
  data2 <- cbind(Type ="knn_loess", test2$PMAD_sd)
  
  test3 <- Total_data_PMAD(knn.rlr.dat)
  data3 <- cbind(Type ="knn_rlr", test3$PMAD_sd)
  
  test4 <- Total_data_PMAD(lls.vsn.dat)
  data4 <- cbind(Type ="lls_vsn", test4$PMAD_sd)
  
  test5 <- Total_data_PMAD(lls.loess.dat)
  data5 <- cbind(Type ="lls_loess", test5$PMAD_sd)
  
  test6 <- Total_data_PMAD(lls.rlr.dat)
  data6 <- cbind(Type ="lls_rlr", test6$PMAD_sd)
  
  test7 <- Total_data_PMAD(svd.vsn.dat)
  data7 <- cbind(Type ="svd_vsn", test7$PMAD_sd)
  
  test8 <- Total_data_PMAD(svd.loess.dat)
  data8 <- cbind(Type ="svd_loess", test8$PMAD_sd)
  
  test9 <- Total_data_PMAD(svd.rlr.dat)
  data9 <- cbind(Type ="svd_rlr", test9$PMAD_sd)
  
  #Combining all the above results
  total_pmad_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_sd <-total_pmad_overall_sd2%>%dplyr::slice_min(PMAD_sd, n=1, with_ties = TRUE)
  original_cols <- c("Overall_Type.PMAD_sd", "Overall_value.PMAD_sd")
  colnames(total_pmad_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PMAD<- cbind(final_Group_data_PMAD_mean, final_Group_data_PMAD_median,
                      final_Group_data_PMAD_sd, total_pmad_overall_mean,
                      total_pmad_overall_median, total_pmad_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PMAD), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PMAD_names1 <- result_PMAD[, -col_indices]
  n <- matrix(t(result_PMAD_names1), ncol=1)
  all_elements <- unlist(n)
  
  # Get the frequency of each element in the dataframe
  freq <- table(all_elements)
  
  # Find the most occurring element
  PMAD_best_combination <- names(freq)[which.max(freq)]
  
  #Finding the best combination
  Best_combinations <- cbind(PCV_best_combination, PEV_best_combination, PMAD_best_combination)
  
  #Adding names to table
  Combinations <- c("knn_vsn", "knn_loess", "knn_rlr",
                    "lls_vsn", "lls_loess", "lls_rlr",
                    "svd_vsn", "svd_loess", "svd_rlr")
  #Extracting PCV values for all combinations
  pcv_group_mean <- total_Group_data_PCV_mean
  pcv_group_median <- total_Group_data_PCV_median
  pcv_group_sd <- total_Group_data_PCV_sd
  
  pcv_overall_mean <- as.data.frame(total_pcv_overall_mean2[,-1])
  colnames(pcv_overall_mean) <- "Overall_PCV_mean"
  pcv_overall_median <- as.data.frame(total_pcv_overall_median2[,-1])
  colnames(pcv_overall_median) <- "Overall_PCV_median"
  pcv_overall_sd <- as.data.frame(total_pcv_overall_sd2[,-1])
  colnames(pcv_overall_sd) <- "Overall_PCV_sd"
  PCV_table2 <- suppressMessages(cbind(pcv_group_mean, pcv_group_median, pcv_group_sd,
                                       pcv_overall_mean, pcv_overall_median, pcv_overall_sd))
  PCV_table1 <- PCV_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PCV_table <- cbind(Combinations, PCV_table1)
  
  #Extracting PEV values for all combinations
  pev_group_mean <- total_Group_data_PEV_mean
  pev_group_median <- total_Group_data_PEV_median
  pev_group_sd <- total_Group_data_PEV_sd
  
  pev_overall_mean <- as.data.frame(total_pev_overall_mean2[,-1])
  colnames(pev_overall_mean) <- "Overall_PEV_mean"
  pev_overall_median <- as.data.frame(total_pev_overall_median2[,-1])
  colnames(pev_overall_median) <- "Overall_PEV_median"
  pev_overall_sd <- as.data.frame(total_pev_overall_sd2[,-1])
  colnames(pev_overall_sd) <- "Overall_PEV_sd"
  PEV_table2 <- suppressMessages(cbind(pev_group_mean, pev_group_median, pev_group_sd,
                                       pev_overall_mean, pev_overall_median, pev_overall_sd))
  PEV_table1 <- PEV_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PEV_table <- cbind(Combinations, PEV_table1)
  
  #Extracting PMAD values for all combinations
  pmad_group_mean <- total_Group_data_PMAD_mean
  pmad_group_median <- total_Group_data_PMAD_median
  pmad_group_sd <- total_Group_data_PMAD_sd
  
  pmad_overall_mean <- as.data.frame(total_pmad_overall_mean2[,-1])
  colnames(pmad_overall_mean) <- "Overall_PMAD_mean"
  pmad_overall_median <- as.data.frame(total_pmad_overall_median2[,-1])
  colnames(pmad_overall_median) <- "Overall_PMAD_median"
  pmad_overall_sd <- as.data.frame(total_pmad_overall_sd2[,-1])
  colnames(pmad_overall_sd) <- "Overall_PMAD_sd"
  PMAD_table2 <- suppressMessages(cbind(pmad_group_mean, pmad_group_median, pmad_group_sd,
                                        pmad_overall_mean, pmad_overall_median, pmad_overall_sd))
  PMAD_table1 <- PMAD_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PMAD_table <- cbind(Combinations, PMAD_table1)
  
  result_list <- list("Best combinations" = as.data.frame(Best_combinations), "PCV Result" = PCV_table, "PEV Result" = PEV_table, "PMAD Result" =PMAD_table,
                      "vsn_data" = cbind(com_data_ID, vsn.dat), "loess_data" = cbind(com_data_ID, loess.dat), "rlr_data" = cbind(com_data_ID, rlr.dat),
                      "knn_vsn_data" = cbind(com_data_ID,knn.vsn.dat),  "knn_loess_data" = cbind(com_data_ID,knn.loess.dat),
                      "knn_rlr_data" = cbind(com_data_ID,knn.rlr.dat), "lls_vsn_data" = cbind(com_data_ID, as.data.frame(lls.vsn.dat)),
                      "lls_loess_data" =  cbind(com_data_ID, as.data.frame(lls.loess.dat)), "lls_rlr_data" = cbind(com_data_ID, as.data.frame(lls.rlr.dat)),
                      "svd_vsn_data" =  cbind(com_data_ID, as.data.frame(svd.vsn.dat)), "svd_loess_data" =  cbind(com_data_ID, as.data.frame(svd.loess.dat)),
                      "svd_rlr_data" =  cbind(com_data_ID, as.data.frame(svd.rlr.dat)))
  sink()
  return (result_list)
}
