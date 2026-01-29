# Functions used for conducting changepoint analysis.
# Various internal and exported fns.

#' Run a Monte-Carlo changepoints analysis on supplied paleoclimate record data
#'
#' @description Apply MCCPT to paleoclimate data
#'
#' @param sites_data formatted paleoclimate data, output from import_files.
#' @param age_lowerbound minimum age for the analysis in years BP.
#' @param age_upperbound maximum age for the analysis in years BP.
#' @param output_folder full file path to a folder where the outputs will be save.
#' @param minseg_len NULL or specify a minimum segment length. User will be prompted for each record if none is specified.
#' @param n_cpts number of changepoints to attempt. User will be prompted on a per-record basis if none is specified.
#' @param cpt_calc which type of changepoint calculation? One of "meanvar" or "mean".
#' @param cpt_method passed to method in cpt.mean or cpt.meanvar. See \link[changepoint]{cpt.mean}
#' @param penalty_value value used for pen.value in cpt.mean or cpt.meanvar.
#' @param penalty_type character; value passed to penalty in cpt.mean or cpt.meanvar.
#' @param save TRUE/FALSE to save outputs to the output_folder.
#' @param rtn TRUE/FALSE to return changepoint data.
#' @param rev_y TRUE/FALSE to reverse y axis on plots
#' @param uncertainty_res Resolution in years for uncertainty resolution within the PrC calculation. Adjust this for shorter series.
#'
#' @importFrom openxlsx addWorksheet
#' @importFrom openxlsx saveWorkbook
#' @importFrom openxlsx writeData
#' @importFrom openxlsx createWorkbook
#'
#' @export
#'
conduct_MCCPT <- function(sites_data,
                          age_lowerbound = 6000,
                          age_upperbound = 22000,
                          output_folder = paste0(getwd(),"/CPT_outputs/"),
                          minseg_len = NULL,
                          n_cpts = NULL,
                          cpt_calc = "meanvar",
                          cpt_method = "BinSeg",
                          penalty_value = "2*log(n)",
                          penalty_type = "Manual",
                          save = TRUE,
                          rtn = TRUE,
                          rev_y = FALSE,
                          verbose = TRUE,
                          uncertainty_res = 20){
  ## Input handling
  # For input data.
  if(is.list(sites_data)){
    if(verbose){
      message("sites_data detected as list. Format will be assumed as correct, i.e. an output from import_files()")
    }
  } else if(is.character(sites_data)){
    sites_data <- import_files(sites_data)
  }
  # For output folder
  if(!dir.exists(output_folder)){
    dir.create(output_folder)
  }
  ## Setup for main loop
  # First, a data frame containing status information for each site
  model_status_df <- data.frame(matrix(ncol = 6, nrow = length(sites_data)))
  colnames(model_status_df) <- c("site_name","model_type","reason_for_failure", "minseg_length", "upperbound", "lowerbound")
  model_status_df[,1] <- names(sites_data)
  ## Init the outputs list
  outputs_list <- vector('list', length(sites_data)) %>%
    'names<-'(c(names(sites_data)))
  ## Start main loop
  # Initialise cpt list
  sites_list <- vector('list',length(sites_data))
  for (i in seq_along(sites_list)){
    # i = 1
    this_site_outputs <- vector('list',length = 1) %>%
      'names<-'(c('prcurves'))
    # message will give iteration number in the console.
    if(verbose){
      message("Running CPT analysis for site ",i,"/",nrow(sites_data)," named ",names(sites_data)[i])
    }
    # Message
    if(verbose & sites_data[[i]]$metadata[which(sites_data[[i]]$metadata[,1] == "Data type"),2] == "Compositional"){
      message("Generating principal curve scores for compositional data types.")
    } else if(verbose & sites_data[[i]]$metadata[which(sites_data[[i]]$metadata[,1] == "Data type"),2] == "Single"){
      message("Analysing as single-series data type.")
    }
    # Run PrC
    PrC_results <- generate_PrC(site_data = sites_data[[i]],
                                age_upper = age_upperbound,
                                age_lower = age_lowerbound)
    dat_i = PrC_results$scrs
    time_i = PrC_results$time
    ageits_i <- PrC_results$ageits
    # model_status_df = status_df
    # message
    if(verbose){
      message("Running changepoint model.")
    }
    cpt.plot <- run_cpts(site_data = sites_data[[i]],
                         site_name = names(sites_data)[i],
                         age_lowerbound = age_lowerbound,
                         age_upperbound = age_upperbound,
                         minseg_len = NULL,
                         n_cpts = NULL,
                         cpt_calc = cpt_calc,
                         cpt_method = cpt_method,
                         penalty_value = penalty_value,
                         penalty_type = penalty_type,
                         PrC_results = PrC_results,
                         status_df = model_status_df,
                         rev_y = rev_y,
                         uncertainty_res = uncertainty_res,
                         verbose = verbose)

    if(verbose){
      message("Extracting age model iteration changepoint densities")
    }
    # This returns a list of density data associated with each changepoint.
    density_list <- CPTwindow(cpt.plot = cpt.plot, ageits_i = ageits_i, plot=FALSE)
    ncpts = length(cpt.plot@cpts)
    if(verbose){
      message("Plotting changepoint densities")
    }
    # This applies over the density list, extracts the y values of each, and finds the max.
    maxy <- max(unlist(lapply(density_list,"[[",'y')), na.rm = T)
    # Loop through and plot changepoint positions on iteration densities
    for(c in seq_along(density_list)){
      if(!is.null(density_list[[c]])){
        par(fig=c(0.5, 1, 0.2, 0.9), mar=c(2,2,2,1), new=TRUE)
        plot(density_list[[c]], type="n", ylab="", xlab="", xlim=c(age_lowerbound, age_upperbound), ylim=c(0,maxy),  main="", cex.axis=1.2)
        polygon(density_list[[c]]$x, density_list[[c]]$y, col="grey", border=NA)
        CPTwindow(cpt.plot = cpt.plot, ageits = ageits_i, plot=FALSE)
        mtext("Age (cal. yr BP)", side=1, cex=1.3, line=3)
        mtext(paste0(names(sites_data)[i]," PDFs"), side=3, adj=0, cex=1.3, line=0.2)
        box(lwd=2)
      }
    }
    ## Save plot
    if(save){
      dev.copy(pdf,paste0(output_folder,names(sites_data)[i],"_CPT.pdf"),  width=10, height=6)
      if (length(dev.list()!=0)) {dev.off()}
      dev.copy(png,paste0(output_folder,names(sites_data)[i],"_CPT.png"),  width=10*100, height=6*100)
      if (length(dev.list()!=0)) {dev.off()}
    }
    ## Output results
    to_keep <- c("x","y","mean","median","wgtmean","firstQ","thirdQ")
    return_list_elements <- function(list, nms){
      elements <- grep(paste(to_keep, collapse = "|"), names(density_list[[1]]))
      list_ret <- list[elements]
    }
    save_list <- vector('list',length = length(density_list))
    for(c in seq_along(density_list)){
      if(!is.null(density_list[[c]])){
        save_list[[c]] <- return_list_elements(density_list[[c]],to_keep)
        save_list[[c]][["thirdQ"]] <- save_list[[c]][["thirdQ"]][["3rd Qu."]]
        save_list[[c]][["firstQ"]] <- save_list[[c]][["firstQ"]][["1st Qu."]]
      }
    }
    # Remove empty final changepoint left over
    save_list <- save_list[lengths(save_list) != 0]
    # Save, if save == TRUE.
    if(save){
      # Get save filename
      savename <- paste0(output_folder,names(sites_data)[i],"-MCCPTs.xlsx")
      # init workbook, create sheets, then write data-frame-coerced output data to the sheets.
      if(verbose){
        message("Saving data to: ", output_folder)
      }
      wb <- createWorkbook()
      for(s in seq_along(save_list)){
        addWorksheet(wb,paste0("CPT-",s))
      }
      for(s in seq_along(save_list)){
        save_df <- convert_to_df(save_list[[s]])
        writeData(wb, paste0("CPT-",s), save_df, startRow = 1, startCol = 1)
      }
      # Save!
      saveWorkbook(wb, file = savename, overwrite = TRUE)
    }
    # If the user has specified rtn as true, put the cpt data into the output list.
    if(rtn){
      outputs_list[[i]] <- save_list
    }
    # Finish loop
    if (length(dev.list()!=0)) {dev.off()}
    message("completed site ",i,"/",length(sites_data),"\n","-------")
  }
  # Return per-site changepoints data
  if(rtn){
    return(outputs_list)
  }
}

#' Run a changepoints analysis (internal)
#'
#' @description
#' Function used to run a changepoint analysis on a given set of data.
#'
#' @param site_data data from one site, a subset of a list returned by import_data()
#' @param site_name the name of the current site.
#' @param age_lowerbound lower age bound
#' @param age_upperbound upper age bound
#' @param minseg_len NULL or a numeric value for minimum number of segments
#' @param n_cpts NULL or a numeric value for the number of changepoints
#' @param cpt_calc which type of changepoint calculation? One of "meanvar" or "mean".
#' @param cpt_method passed to method in cpt.mean or cpt.meanvar. See \link[changepoint]{cpt.mean}
#' @param penalty_value value used for pen.value in cpt.mean or cpt.meanvar.
#' @param penalty_type character; value passed to penalty in cpt.mean or cpt.meanvar.
#' @param PrC_results results from generate_PrC()
#' @param rev_y TRUE/FALSE to reverse y axis
#' @param uncertainty_res resolution of uncertainty in years
#' @param verbose TRUE/FALSE to print status as fn proceeds
#'
#' @importFrom changepoint cpt.meanvar
#'
#' @noRd
#'
run_cpts <- function(site_data,
                     site_name,
                     age_lowerbound,
                     age_upperbound,
                     minseg_len = NULL,
                     n_cpts = NULL,
                     cpt_calc = "meanvar",
                     cpt_method = "BinSeg",
                     penalty_value = "2*log(n)",
                     penalty_type = "Manual",
                     PrC_results = NULL,
                     status_df = NULL,
                     rev_y = NULL,
                     uncertainty_res = NULL,
                     verbose){
  # Data passover
  dat_i = PrC_results$scrs
  time_i = PrC_results$time
  model_status_df = status_df
  # Get n observation
  obvs<-as.numeric(nrow(as.data.frame(time_i)))
  # Here k opeates
  likemodel=0 # set k to zero so that you can reset the loop
  while(likemodel<1){
    # Generate initial plot
    plot_cpt_pre(PrC_results,
                 age_upperbound = age_upperbound,
                 age_lowerbound = age_lowerbound,
                 rev_y = rev_y,
                 uncertainty_res = uncertainty_res,
                 verbose = verbose,
                 name = site_name)
    ## Segment length selection
    if(is.null(minseg_len)){
      message("Set MinSeg length for ",site_name,
              "\n","Number of observations between ",age_lowerbound,"-",age_upperbound," is ", obvs)
      k <- readline(prompt = "MinSeg = ")
      k<-as.integer(k)
    } else {
      if(!is.numeric(minseg_len)){stop("minseg_len must be NULL or numeric")}
      k = minseg_len
    }
    ## NB: model status currently not working
    # model_status_df[i,4] <- k
    ## N changepoint selection
    if(is.null(n_cpts)){
      # Q = max n changepoints to search for
      message("Set n changepoints for ", site_name)
      q <- readline(prompt = "Q = ")
      q<-as.integer(q)
    } else {
      if(!is.numeric(n_cpts)){stop("n_cpts must be NULL or numeric")}
      q <- n_cpts
    }
    # first we set the skip logical to FALSE. If it fails the trycatch,
    # i.e. the model is so bad it fails, then the skip to next is set to true
    # and the user is prompted to try a different minseg and Q.
    skip_to_next <- FALSE
    # Attempt model, output error to a function (no real idea why R does this)
    if(cpt_calc == "meanvar"){
      test <- tryCatch(cpt.meanvar(dat_i, penalty = penalty_type, pen.value = penalty_value,
                                   minseglen = k, Q = q, method = cpt_method, class = TRUE),
                       error = function(e) {
                         skip_to_next <<- TRUE
                       })
    } else if (cpt_calc == "mean"){
      test <- tryCatch(cpt.mean(dat_i, penalty = penalty_type, pen.value = penalty_value,
                                   minseglen = k, Q = q, method = cpt_method, class = TRUE),
                       error = function(e) {
                         skip_to_next <<- TRUE
                       })
    }
    # if skip to next is true, skip iteration.
    if(skip_to_next) {
      message("\n","Models failed due to minseglen error. Skipping iteration.","\n","------")
      # model_status_df[i,2] = "FAIL"
      # model_status_df[i,3] = "minseglenth error"
      midY <- median(c(min(dat_i), max(dat_i)))
      lines(c(age_lowerbound, age_upperbound), c(midY,midY), col="blue", lwd=2)
      mtext(paste0(site_name), side=3, adj=0, cex=1.3, line=0.2)
      mtext("Age (cal. yr BP)", side=1, cex=1.3, line=3)
      box(lwd=2)
      # dev.copy(pdf,paste0(output_folder,site_name"_CPT.pdf"),  width=10, height=6)
      # if (length(dev.list()!=0)) {dev.off()}
      # message("completed site ",i,"/",nrow(sites_df),"\n","-------")
      next
    } else {
      if(cpt_calc == "meanvar"){
        bm1_i <- cpt.meanvar(dat_i, penalty = penalty_type, pen.value = penalty_value,
                             minseglen = k, Q = q, method = cpt_method, class = TRUE )
      } else if(cpt_calc == "mean"){
        bm1_i <- cpt.mean(dat_i, penalty = penalty_type, pen.value = penalty_value,
                             minseglen = k, Q = q, method = cpt_method, class = TRUE )
      }
      # bm1_i <- cpt.meanvar(dat_i, penalty = "Manual",pen.value="2*log(n)",
      #                      minseglen = k, Q=q, method="BinSeg", class = TRUE )
      ncpts<-as.numeric(nrow(as.data.frame(bm1_i@cpts)))
      cpt.plot<-bm1_i
      plot.cpts(cpt.plot, time=time_i, timescale = "BP")
      likemodel <- readline("Do you like this model: 0=No, 1=Yes ")
      if(likemodel == 0){

      }
      #k=1
    }
  }
  if (skip_to_next){
    next
  }
  cpt.plot
}

#' Generate principle curve scores for a given set of data
#'
#' @description
#' Do the PrC
#'
#' @param site_data data from a single
#' @param age_upper upperbound of age
#' @param age_lower lowerbound of age
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom analogue prcurve
#' @importFrom vegan scores
#'
#' @noRd
#'
generate_PrC <- function(site_data,
                         age_upper = age_upperbound,
                         age_lower = age_lowerbound,
                        output_folder = paste0(getwd(),"/CPT_outputs/"),
                          save = TRUE){
  ## var passover
  age_upperbound = age_upper
  age_lowerbound = age_lower
  ## Age handling
  # See e.g. data for format.
  bestage_i <- site_data$ages$Mean
  ages_i <- site_data$ages %>%
    select(-c(Median)) %>%
    select(-c(which(grepl('Depth', colnames(.)))))
  ## Per-iteration age handling
  if(is.null(age_lowerbound) & is.null(age_upperbound)){
    message("Set age bounds for site ",i,"/",nrow(sites_df)," named ",sites_df[i,1])#,
    # "\n","Number of observations is ", obvs)
    # Prompt user for age lowerbound
    l <- readline(prompt = "Lowerbound = ")
    age_lowerbound_i <- as.integer(l)
    model_status_df$upperbound[i] <- age_lowerbound_i
    # Do the same for upperbound
    u <- readline(prompt = "Uowerbound = ")
    age_upperbound_i <- as.integer(u)
    model_status_df$upperbound[i] <- age_upperbound_i
  } else {
    # This is the default; user should probably supply age bounds to be consistent.
    age_lowerbound_i = age_lowerbound
    age_upperbound_i = age_upperbound
  }
  ## Pollen site iteration
  PrC_res = vector('list', length = 3) %>%
    'names<-'(c('scrs','time','ageits'))
  # Compositional data only (typically pollen or diatoms)
  if(site_data$metadata[which(site_data$metadata$category == 'Data type'),2] == "Compositional" ){
    # Reset the graphics device
    if (length(dev.list()!=0)) {dev.off()}
    # extract pollen % data from imported file
    dat_i = site_data$data [,-c(1)] %>%
      mutate(age = bestage_i) %>%
      relocate(age)
    # subset data to desired age range, specified as input parameters age_lowerbound and age_upperbound.
    dat_i <- subset(dat_i, dat_i$age > age_lowerbound_i & dat_i$age < age_upperbound_i)
    # This line would be the one to tweak if you wanted to ensure data points outside of the 15-35k range are included
    ages_i <- subset(ages_i, ages_i$Mean > age_lowerbound_i & ages_i$Mean < age_upperbound_i)
    bestage_i <- subset(bestage_i, bestage_i > age_lowerbound_i & bestage_i < age_upperbound_i)
    # Delete final column (contains 'Rumex' data)
    pollen_i <- dat_i %>%
      select(-c(age))
    # Run a prcurve on the pollen data
    pollen_pc_i <- prcurve(
      pollen_i,
      method = 'pca',
      trace = TRUE,
      plotit = FALSE,
      vary = FALSE,
      penalty = 1.2)
    # system sleep so that the output shows up.
    # Sys.sleep(0)
    # store prcurve scores
    scrs_i <- scores(pollen_pc_i)
    # Sys.sleep(0)
    PrC_res$scrs <- as.numeric(na.omit(scrs_i))
    PrC_res$time <- bestage_i
    PrC_res$ageits <- ages_i
    # dat_i <- as.numeric(na.omit(scrs_i))
    # time_i <- bestage_i
    # ageits_i <- ages_i
  } else {
    # NON-POLLEN SITES
    # Reset null device
    if (length(dev.list()!=0)) {dev.off()}
    # Get dat
    dat_i = site_data$data %>% select(-c(1)) %>%
      mutate(age = bestage_i) %>%
      relocate(age)
    # subset data to desired age range, specified as input parameters age_lowerbound and age_upperbound.
    dat_i <- subset(dat_i, dat_i$age > age_lowerbound_i & dat_i$age < age_upperbound_i)
    # This line would be the one to tweak if you wanted to ensure data points outside of the 15-35k range are included
    ages_i <- subset(ages_i, ages_i$Mean > age_lowerbound_i & ages_i$Mean < age_upperbound_i)
    bestage_i <- subset(bestage_i, bestage_i > age_lowerbound_i & bestage_i < age_upperbound_i)
    # dat_i <- read.xlsx(paste0(data_dir,"/",sites_df[i,1],".xlsx"), sheetName = "Data", header = TRUE)
    # replace sample with corresponding age
    dat_i[,1] <- bestage_i
    colnames(dat_i)[1] <- "age"
    # subset data to desired age range, specified as input parameters age_lowerbound and age_upperbound.
    dat_i <- subset(dat_i, dat_i$age > age_lowerbound_i & dat_i$age < age_upperbound_i)
    # store prcurve scores
    scrs_i <- as.matrix(dat_i[,2])
    PrC_res$scrs <- as.numeric(na.omit(scrs_i))
    PrC_res$time <- bestage_i
    PrC_res$ageits <- ages_i
    # dat_i <- as.numeric(na.omit(scrs_i))
    # time_i <- bestage_i
    # ageits_i <- ages_i
  }
  
  # --- ADDED: save PrC_res if requested (no other changes to logic) ----------------
  if (isTRUE(save)) {
    # Ensure output folder exists
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
    }
    # Build a simple, robust filename
    # Try to pull a site-like name from metadata; fall back to "site"
    site_stub <- tryCatch({
      md <- site_data$metadata
      if (!is.null(md) && is.data.frame(md) && nrow(md) > 0) {
        # Look for likely name fields; adjust if you have a known key
        idx <- which(tolower(md$category) %in% c("site", "site name", "name", "id"))
        if (length(idx) > 0) as.character(md[idx[1], 2]) else "site"
      } else {
        "site"
      }
    }, error = function(e) "site")
    site_stub <- gsub("[^[:alnum:]_\\-]+", "_", site_stub)

    # Include bounds in filename; use the already computed _i variables
    rds_path <- file.path(
      output_folder,
      paste0(site_stub, "_PrC_", age_lowerbound_i, "-", age_upperbound_i, "BP.rds")
    )

    saveRDS(PrC_res, rds_path)
    message("Saved PrC_res to: ", rds_path)
  }
  # --- END ADDED ------------------------------------------------------------------

  # Return results
  return(PrC_res)
}

