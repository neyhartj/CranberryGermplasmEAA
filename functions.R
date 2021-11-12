## Project name here
##
## This is a script with ad hoc functions for this project
##

# # Read the genomic data into sambada
# fileName = file.path(data_dir, "gcWild_filtered_phased_markers.vcf")
# outputFile = file.path(data_dir, "gcWild_genos_prepared.csv")
# saveGDS = TRUE
# mafThresh = 0
# missingnessThresh = 1
# interactiveChecks = FALSE
# ldThresh = NULL
# mgfThresh = NULL
# directory = NULL
# interactiveChecks = FALSE
# verbose = FALSE


prepareGeno2 <- function (fileName, outputFile, saveGDS, mafThresh = NULL, missingnessThresh = NULL,
                          ldThresh = NULL, mgfThresh = NULL, directory = NULL, interactiveChecks = FALSE,
                          verbose = FALSE) {

  requireNamespace("R.SamBada")

  if (typeof(fileName) != "character")
    stop("fileName argument supposed to be character")
  if (!file.exists(fileName))
    stop("Input file not found.")
  if (typeof(outputFile) != "character")
    stop("outputFile argument supposed to be character")
  extensionO = substr(outputFile, gregexpr("\\.", outputFile)[[1]][length(gregexpr("\\.",
                                                                                   outputFile)[[1]])] + 1, nchar(outputFile))
  if (extensionO != "csv")
    stop("outputFile must have a .csv extension")
  if (typeof(saveGDS) != "logical")
    stop("saveGDS argument supposed to be logical")
  if (!is.null(mafThresh)) {
    if (typeof(mafThresh) != "double")
      stop("mafThresh argument supposed to be decimal number")
    if (mafThresh > 1 | mafThresh < 0)
      stop("mafThresh argument supposed to be between 0 and 1")
  }
  if (!is.null(missingnessThresh)) {
    if (typeof(missingnessThresh) != "double")
      stop("missingnessThresh argument supposed to be decimal number")
    if (missingnessThresh > 1 | missingnessThresh < 0)
      stop("missingnessThresh argument supposed to be between 0 and 1")
  }
  if (!is.null(mgfThresh)) {
    if (typeof(mgfThresh) != "double")
      stop("mgfThresh argument supposed to be decimal number")
    if (mgfThresh > 1 | mgfThresh < 0)
      stop("mgfThresh argument supposed to be between 0 and 1")
  }
  if (!is.null(ldThresh)) {
    if (typeof(ldThresh) != "double")
      stop("ldThresh argument supposed to be decimal number")
    if (ldThresh > 1 | ldThresh < 0)
      stop("ldThresh argument supposed to be between 0 and 1")
  }
  if (typeof(interactiveChecks) != "logical")
    stop("interactiveChecks argument supposed to be logical")
  if (typeof(verbose) != "logical")
    stop("verbose argument supposed to be logical")
  fileName_base = basename(fileName)
  filename_short = substr(fileName, 1, gregexpr("\\.",
                                                fileName)[[1]][length(gregexpr("\\.", fileName)[[1]])] -
                            1)
  filename_short2 = substr(fileName_base, 1, gregexpr("\\.",
                                                      fileName_base)[[1]][length(gregexpr("\\.", fileName_base)[[1]])] -
                             1)
  extension = substr(fileName, gregexpr("\\.", fileName)[[1]][length(gregexpr("\\.",
                                                                              fileName)[[1]])] + 1, nchar(fileName))
  if (!is.null(directory)) {
    changePath(directory)
  }
  tryCatch(suppressWarnings(system2("recode-plink", wait = TRUE,
                                    stdout = FALSE, stderr = FALSE)), error = function(e) {
                                      stop("sambada's recode-plink is not available. You should first download sambada and either put the binary folder to the path environmental variable or specify the path in the directory input argument")
                                    })
  if (extension == "ped" & is.null(mafThresh) & is.null(missingnessThresh) &
      is.null(mgfThresh) & is.null(ldThresh)) {
    if (!file.exists(paste(filename_short, ".map",
                           sep = "")))
      stop(".map input file not found. Same name as .ped mandatory")
    numSNP = nrow(read.table(paste(filename_short, ".map",
                                   sep = ""), colClasses = c("character",
                                                             rep("NULL", 3)), sep = "\t"))
    f = file(fileName, open = "rb")
    numIndiv = 0L
    while (length(chunk <- readBin(f, "raw", 20000)) >
           0) {
      numIndiv = numIndiv + sum(chunk == as.raw(10L))
    }
    close(f)
    system2("recode-plink", args = c(numIndiv, numSNP,
                                     filename_short, outputFile))
    return(NA)
  }
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("Package \"SNPRelate\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Package \"gdsfmt\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  tmp = tempdir()
  current = getwd()
  if (verbose == TRUE) {
    print("Creating and opening GDS file")
    print("===================================")
  }
  if (saveGDS == FALSE) {
    outputDir = tmp
  } else {
    outputDir = getwd()
  }
  if (interactiveChecks == TRUE & extension != "gds" &
      file.exists(paste0(filename_short, ".gds"))) {
    useGDS = readline(prompt = "A .gds file has been found with the same name. Would you like to use it (press y) or rewrite a new file (press any other key): ")
    if (useGDS == "y") {
      gds_file = paste0(filename_short, ".gds")
    }
    else {
      gds_file = createGDSfile(fileName, outputDir)
    }
  } else {
    gds_file = R.SamBada:::createGDSfile(fileName, outputDir)
  }
  gds_obj = SNPRelate::snpgdsOpen(gds_file)
  on.exit(SNPRelate::snpgdsClose(gds_obj))
  MGF <- NULL
  if (is.null(mgfThresh) == FALSE) {
    Rcpp::cppFunction("\n      NumericVector MGF(NumericMatrix xx, double maxMGFAllowed){\n        int xrow = xx.nrow() ;\n        int xcol = xx.ncol();\n        int aa;\n        int Aa;\n        int AA;\n        int sum_max;\n        int mgf;\n        NumericVector yy(xcol);\n        NumericVector yybool(xcol);\n        NumericVector yysnpid(xcol);\n        int k=0;\n        for(int i = 0; i < xcol; i++) {\n        aa=0;\n        Aa=0;\n        AA=0;\n        for(int j = 0; j < xrow; j++){\n        if(xx(j,i)==0){\n        aa++;\n        }\n        else if(xx(j,i)==1){\n        Aa++;\n        }\n        else if(xx(j,i)==2){\n        AA++;\n        }\n        }\n        if(aa>=Aa){\n        sum_max=aa;\n        if(AA>aa){\n        sum_max=AA;\n        }\n        }\n        else{\n        if(AA>=Aa){\n        sum_max=AA;\n        }\n        else{\n        sum_max=Aa;\n        }\n        }\n        if(aa+AA+Aa>0){\n        mgf=(sum_max*100)/(aa+Aa+AA);\n        }\n        else{\n        mgf=0;\n        }\n        yy(i)=mgf;\n        if(mgf>maxMGFAllowed*100){\n        yybool(i)=0;\n        \n        }\n        else{\n        yybool(i)=1;\n        yysnpid(k)=i+1;\n        k++;\n        }\n        \n        }\n        \n        //NumericVector yysnpid2(k);\n        //yysnpid2(yysnpid.begin() , yysnpid.begin() + k);\n        if(maxMGFAllowed>=0){\n        return yysnpid;\n        }\n        else{\n        return yy;\n        }\n    }")
  }
  if (verbose == TRUE) {
    print("Filtering using SNPRelate in process")
    print("===================================")
  }
  if (interactiveChecks == TRUE) {
    snpSummary = SNPRelate::snpgdsSNPRateFreq(gds_obj, with.snp.id = TRUE)
    if (is.null(mafThresh) == FALSE) {
      hist(snpSummary$MinorFreq, breaks = 100, xlab = "Minor allele Frequency",
           main = "Histogram of minor allele frequency",
           xlim = c(0, max(mafThresh, max(snpSummary$MinorFreq))))
      abline(v = mafThresh, col = "red")
      mafThresh2 = readline(prompt = "Would you like to change your maf threshold? (press n if no, or enter a new threshold): ")
      if (grepl("[[:digit:]\\.-]", mafThresh2)) {
        if (as.numeric(mafThresh2) > 1 | as.numeric(mafThresh2) <
            0)
          stop("mafThresh argument supposed to be between 0 and 1")
        mafThresh = as.numeric(mafThresh2)
      }
    }
    if (is.null(missingnessThresh) == FALSE) {
      hist(snpSummary$MissingRate, breaks = 100, xlab = "Missingness",
           main = "Histogram of missingness", xlim = c(0,
                                                       max(missingnessThresh, max(snpSummary$MissingRate))))
      abline(v = missingnessThresh, col = "red")
      missingnessThresh2 = readline(prompt = "Would you like to change your missingness threshold? (press n if no, or enter a new threshold): ")
      if (grepl("[[:digit:]\\.-]", missingnessThresh2)) {
        if (as.numeric(missingnessThresh2) > 1 | as.numeric(missingnessThresh2) <
            0)
          stop("missingnessThresh argument supposed to be between 0 and 1")
        missingnessThresh = as.numeric(missingnessThresh2)
      }
    }
    if (is.null(ldThresh) == FALSE) {
      if (length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj,
                                                      "snp.id"))) > 1e+06) {
        compute.LD = readline(prompt = "Would you like to compute the graph of discarded SNP vs LD (takes a long time given the number of SNPs) ? (press y for yes, any other key for no): ")
      }
      else {
        compute.LD = "y"
      }
      if (compute.LD == "y") {
        snp_pruned = vector()
        i = 1
        for (ld in seq(0, 1, 0.05)) {
          snp_pruned[i] = length(unlist(SNPRelate::snpgdsLDpruning(gds_obj,
                                                                   ld.threshold = ld, verbose = FALSE, remove.monosnp = FALSE)))
          i = i + 1
        }
        barplot(snp_pruned, names.arg = seq(0, 1, 0.05),
                space = 0, col = "white", xlab = "LD",
                ylab = "Number of SNPs pruned", main = "Histogram of number of pruned SNP by LD")
        abline(v = ldThresh * length(seq(0, 1, 0.05)) +
                 1 - ldThresh, col = "red")
        ldThresh2 = readline(prompt = "Would you like to change your LD threshold? (press n if no, or enter a new threshold): ")
        if (grepl("[[:digit:]\\.-]", ldThresh2)) {
          if (as.numeric(ldThresh2) > 1 | as.numeric(ldThresh2) <
              0)
            stop("ldThresh argument supposed to be between 0 and 1")
          ldThresh = as.numeric(ldThresh2)
        }
      }
    }
    if (is.null(mgfThresh) == FALSE) {
      geno = SNPRelate::snpgdsGetGeno(gds_obj)
      mgf_freq = MGF(geno, -1)
      hist(mgf_freq/100, breaks = 100, main = "Histogram of MGF")
      abline(v = mgfThresh, col = "red")
      mgfThresh2 = readline(prompt = "Would you like to change your MGF threshold? (press n if no, or enter a new threshold): ")
      if (grepl("[[:digit:]\\.-]", mgfThresh2)) {
        if (as.numeric(mgfThresh2) > 1 | as.numeric(mgfThresh2) <
            0)
          stop("mgfThresh argument supposed to be between 0 and 1")
        mgfThresh = as.numeric(mgfThresh2)
      }
    }
  }
  if (is.null(mafThresh)) {
    mafThresh = NaN
  }
  if (is.null(missingnessThresh)) {
    missingnessThresh = NaN
  }
  if (is.null(ldThresh)) {
    snp_filtered = SNPRelate::snpgdsSelectSNP(gds_obj, autosome.only = FALSE,
                                              maf = mafThresh, missing.rate = missingnessThresh,
                                              verbose = FALSE)
  } else {
    gds_pruned = SNPRelate::snpgdsLDpruning(gds_obj, maf = mafThresh,
                                            missing.rate = missingnessThresh, ld.threshold = ldThresh,
                                            verbose = FALSE)
    snp_filtered = unlist(gds_pruned)
  }
  list_snptype = SNPRelate::snpgdsSNPList(gds_obj)
  snp_filtered2 = list_snptype$snp.id[!(list_snptype$allele %in%
                                          c("A/C", "C/A", "A/G", "G/A",
                                            "A/T", "T/A", "C/G", "G/C",
                                            "C/T", "T/C", "G/T", "T/G",
                                            "A/0", "C/0", "G/0", "T/0"))]
  snp_filtered = snp_filtered[!(snp_filtered %in% snp_filtered2)]
  if (!is.null(mgfThresh)) {
    if (interactiveChecks == FALSE) {
      geno = SNPRelate::snpgdsGetGeno(gds_obj)
    }
    snpMGF = MGF(geno, mgfThresh)
    snpMGF = snpMGF[snpMGF > 0]
    list_snp = snpMGF[snpMGF %in% snp_filtered]
  } else {
    list_snp = snp_filtered
  }
  if (interactiveChecks == TRUE & file.exists(outputFile)) {
    cont = readline(prompt = paste0(outputFile, " already exists and will be replaced. Do you want to continue? (press x to exit, any other key to continue): "))
    if (cont == "x") {
            print("Function ended on user input")
            return(NA)
        }
    }
    SNPRelate::snpgdsGDS2PED(gds_obj, ped.fn = file.path(tmp,
        paste(filename_short2, "_filtered", sep = "")),
        snp.id = list_snp, verbose = FALSE)
    if (verbose == TRUE) {
        print(paste0("Filtering finished: ", length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj,
            "snp.id"))) - length(list_snp), " deleted out of ",
            length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj,
                "snp.id"))), " SNPs"))
        print("=====================================")
    }
    if (verbose == TRUE) {
        print("Running Recodeplink to comply with sambada input file")
        print("=====================================")
    }
    numSNP = length(list_snp)
    numIndiv = length(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds_obj,
        "sample.id")))
    system2("recode-plink", args = c(numIndiv, numSNP,
        file.path(tmp, paste(filename_short2, "_filtered",
            sep = "")), outputFile))
}



# Function to update a formula with a list of predictors
update.formula2 <- function(old, predictors) {
  new <- as.formula(paste0("~ . + ", paste0(predictors, collapse = "+")))
  update.formula(old = old, new = new)
}


# Function to return uniform quantiles given a set of p-values
qqscore <- function(scores) {
  remove <- which(scores == 0)
  scores1 <- if (length(remove) > 0) scores[-remove] else scores

  # Orders
  scores_order <- cbind(scores1, original_order = seq_along(scores1),
                        order = order(scores1, decreasing = TRUE))
  # sort
  scores_order1 <- scores_order[scores_order[,"order"],]

  # Get uniform pvalues
  unif_p <- -log10(ppoints(n = length(scores1)))

  # Return those p-values in the correct order
  unif_p[order(scores_order1[,"original_order"])]

}
