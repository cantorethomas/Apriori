#!/usr/bin/env Rscript

library("optparse")
library("stringr") 

# Implementation of APriori algorithm for pattern mining in databases 
## Finds frequent itemsets using an iterative level-wise approach
## based on candidate generation and evaluation

# Option list ------------------------------------------------------------------

get_params <- function(){
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default=NA,
            help="Dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="patterns.txt",
            help="Output file name [default= %default]", metavar="character"),
    make_option(c("--minsup"), type = "numeric", default=1,
            help="Minimum support for frequent itemsets"),
    make_option(c("--info"), type = "character", default="n",
            help="y/n, print session info"),
    make_option(c("--one"), type = "character", default="n", 
            help="y/n, if just 1-itemsets list has to be returned")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  opt
}


# Check params -----------------------------------------------------------------

check_params <- function(params){
  msg <- c()
  if((!is.na(params[[1]])) | (!file.exists(params[[1]]))) {
    stop("Input file not existing.")
  }
  if(!file.exists(params[[2]])) {
    msg <- c(msg, "No output file given, new one will be created.")
  }
  if(params[[3]] < 0) {
    stop("Negative minsup values are not allowed.")
  }
  if((is.integer(params[[3]])) & (params[[3]] <= 1)) {
    stop("Required 0 < minsup < 1!")
  }
  if (!any(params[[4]] == c("y", "n"))) {
    msg <- c(msg, "Invalid value for --info, default will be used instead.")
  }
  if (params[[4]] == "y") {
    sessionInfo()
  }
  if (!any(params[[5]] == c("y","n"))) {
    msg <- c(msg, "Invalid value for --one, default will be used instead.")
  }
  if (length(msg) != 0) msg
}

# Find frequent 1-itemsets -----------------------------------------------------

find_frequent_1 <- function(data, minsup) { # assuming data as list of itemsets
  freq_vect <- table(unlist(data))
  freq_vect <- freq_vect[freq_vect >= minsup] # min sup as count
  
  if(length(freq_vect) == 0) {
    stop(paste0("No frequent itemsets for --minsup = ", minsup))
  }
  freq_list <- split(names(freq_vect), f = factor(1:length(freq_vect)))
  names(freq_list) <- NULL
  list(freq_list, freq_vect)
}

# has_freq() -------------------------------------------------------------------

has_freq <- function(tmp, freq_list) {
  any(sapply(freq_list, function(x) all((tmp == x))))
}

# apriori_gen() ----------------------------------------------------------------

apriori_gen <- function(freq_list, k) {
  joined_list <- sapply(freq_list, function(x) {
    sapply(freq_list, function(y){
      if(length(x) == 1) { # for 1-itemsets
        tmp <- unique(c(x, y))
        if(length(tmp) != k-1) tmp 
      } else {
        if(all(x[-(k-1)] == y[-(k-1)])){ 
          tmp <- sort(unique(c(x, y)))
          if(length(tmp) != k-1) {
            if(has_freq(tmp[-1], freq_list)) tmp
          }
        }
      }
    }, simplify = T)
  }, simplify = T)
  joined_list <- joined_list[lower.tri(joined_list)]
  joined_list <- joined_list[!sapply(joined_list, is.null)]
  joined_list
}

# get_final() ------------------------------------------------------------------

get_final <- function(cand, freq_tmp){
  d_f <- data.frame(x = unlist(freq_tmp), y = I(cand))
  d_f
}
# Main -------------------------------------------------------------------------

main <- function(data, minsup){
  
  freq_list_1 <- find_frequent_1(data, minsup)
  freq_list <- freq_list_1[[1]]
  
  final_list <- data.frame(freq_list_1[[2]], names(freq_list_1[[2]]))
  if (ncol(final_list) == 3) final_list[1] <- NULL
  colnames(final_list) <- c("x","y")
  
  k <- 2
  cand <- apriori_gen(freq_list, k)
  
  while (length(cand) != 0) {
    freq_tmp <- list() # temporary list for freq collection
    length(freq_tmp) <- length(cand)
    
    names(cand) <- names(freq_tmp) <- 1:length(cand)
    
    for (x in data) { # data scan for frequencies
      for (y in names(cand)) {
        if(all((cand[[y]])%in%(x))) {
          if(is.null(freq_tmp[[y]])) {
            freq_tmp[[y]] <- 1
          } else { 
            freq_tmp[[y]] <- freq_tmp[[y]] + 1
          }
        }
      }
    }
    cand <- cand[!sapply(freq_tmp, is.null)] 
    cand <- cand[unlist(freq_tmp) >= minsup]
    freq_tmp <- freq_tmp[!sapply(freq_tmp, is.null)]
    freq_tmp <- freq_tmp[unlist(freq_tmp) >= minsup]
    names(cand) <- NULL
    final_list <- rbind(get_final(cand, freq_tmp), final_list) # constructing df
    k <- k+1 
    cand <- apriori_gen(cand, k) # new candidates
  }
  final_list
}

# format_out() -----------------------------------------------------------------

format_out <- function(out_df) {
  freqs <- vapply(out_df[,1], function(x) paste0(x,":"), "1:")
  sets <- vapply(out_df[,2], function(x) {
    single_set <- paste(x, collapse = ";")
    single_set
  }, "abc")
  out_df <- data.frame(freqs, sets)
  out_df
}

# Calling ----------------------------------------------------------------------

start <- Sys.time()

# Call for input parameters
params <- get_params()
check_params(params)

# Read in the data
in_df <- scan(params[[1]], what="", sep="\n")
in_df <- strsplit(in_df, ";")

# setting minsup as count 
minsup <- round(params[[3]]*length(in_df[]))

# calling the algorithm 
if (params[[5]] == "y") {
  # creating out_df just for 1-itemsets 
  freq_list_1 <- find_frequent_1(in_df, minsup)
  freq_list <- freq_list_1[[1]]
  out_df <- data.frame(freq_list_1[[2]], names(freq_list_1[[2]]))
  if (ncol(out_df) == 3) out_df[1] <- NULL
} else {
  # running complete algo 
  out_df <- main(in_df, minsup)
}
 

# formatting output 
out_df <- format_out(out_df)

# writing output 
write.table(out_df, 
            params[[2]],
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = "")
end <- Sys.time()

# Out pasting 
paste("Time elapsed: ", round(end-start, 3), " seconds.")
paste("Number of sets found: ", length(out_df[,1]))
paste("Minimum support count: ", minsup)
