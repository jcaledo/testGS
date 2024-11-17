#  install.packages("EnvNJ", type = "source", lib = "/Users/jcaledo/R_pkgs/customLib")

## ---------------------------------------------- ##
#             Testing msa                          #
## ---------------------------------------------- ##
test_that("msas() works properly", {

  skip_on_cran()
  skip_on_travis()

  sequences <- c("MQPIPDVNQRIARISAHLHPPKSQMEESSALRRANCRAKGGAPGFKV",
                 "MSEPIRVLVTGAAGQIAYSLLYSIGNGSVFGKDQPIILVLLDITPMM",
                 "MSEPIRVLVTGAAGQIAYSLLYSIGNGSVFGKDQPIILVLLDITPMM")
  names(sequences) <- c("P19446", "P40925", "P40926")

  a <- msa(sequences = sequences, seqtype = "protein", method = "Muscle")
  b <- msa(sequences = sequences, seqtype = "protein", method = "ClustalOmega")
  c <- msa(sequences = sequences, seqtype = "protein", method = "ClustalW")

  expect_is(a, "data.frame")
  expect_equal(dim(a), c(3,57))
  expect_is(a[[1]], 'character')
  expect_equal(rownames(a), c("P19446", "P40925", "P40926"))

  expect_is(b, "data.frame")
  expect_equal(dim(b), c(3,74))
  expect_is(b[[1]], 'character')
  expect_equal(rownames(b), c("P19446", "P40925", "P40926"))

  expect_is(c, "data.frame")
  expect_equal(dim(c), c(3,49))
  expect_is(c[[1]], 'character')
  expect_true(sum(rownames(c) %in% c("P19446", "P40925", "P40926")) == 3)

})

## ---------------------------------------------------------------- ##
#      msax <- function(sequences, ids, seqtype, method, sfile)       #
## ---------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein, DNA sequences using the msa package.
#' @usage msas(sequences, ids = names(sequences), seqtype = "protein", method = "Muscle", sfile = FALSE)
#' @param sequences vector containing the sequences as strings.
#' @param ids character vector containing the sequences' ids.
#' @param seqtype it should be either "protein" or "dna.
#' @param method currently, "Muscle", ClustalW", "ClustalOmega" are supported.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @details The msa package works without requiring you to install MUSCLE or other alignment tools. The msa package includes its own precompiled binaries for MUSCLE, ClustalW, and Clustal Omega, which are used internally
#' @return Returns an alignment in matrix (dataframe) format. Optionally, the fasta format of the alignment can be saved.
#' @examples \dontrun{msax(sequences = c("APGW", "AGWC", "CWGA"),
#'                        ids = c("a", "b", "c"))}
#' @importFrom msa msa
#' @export

# msax <- function(sequences, ids = names(sequences), seqtype = "protein", method = "Muscle", sfile = FALSE){
#
#   names(sequences) <- ids
#
#   tr <- function(seq, genetic_code = 1){
#     seq <- gsub(" ", "", seq)
#     seq <- strsplit(seq, "")[[1]]
#     output <- paste(translate(seq, numcode = genetic_code), collapse = "")
#
#     return(output)
#   }
#
#   if (length(sequences) < 2) {
#     stop("At least two sequences are required!")
#   } else if (length(sequences) != length(ids)) {
#     stop("The number of sequences and sequences' ids doesn't match!")
#   }
#   if (seqtype == "dna"){
#     dnaSeq <- sequences
#     cod <- strsplit(gsub("(.{3})", "\\1 ", dnaSeq), split = " ")
#     sequences <- unlist(lapply(sequences, function(x) tr(x)))
#   }
#
#   aln <- msa::msa(inputSeqs = sequences, type = seqtype, method = method)
#
#   alnseq <- as.character(aln)
#
#   if (sfile != FALSE){
#     for (i in 1:length(alnseq)){
#       cat(">", ids[i], "\n", alnseq[i], "\n", file = sfile, append = TRUE)
#     }
#   }
#
#   ali <- as.data.frame(do.call(rbind, strsplit(alnseq, split = "")))
#
#   return(ali)
# }





## ---------------------------------------------------------------- ##
#         msa3 <- function(sequences, ids, seqtype, sfile)           #
## ---------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein, DNA or CDS sequences using the R package 'muscle'.
#' @usage msa3(sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE)
#' @param sequences vector containing the sequences as strings.
#' @param seqtype it should be either "prot" of "dna" or "cds" (see details).
#' @param ids character vector containing the sequences' ids.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @details If seqtype is set to "cds" the sequences must not contain stop codons and they will be translated using the standard code. Afterward, the amino acid alignment will be used to lead the codon alignment.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($id) returns the identifiers, the third element ($aln) provides the alignment in fasta format and the fourth element ($ali) gives the alignment in matrix format.
#' @examples \dontrun{msa2(sequences = c("APGW", "AGWC", "CWGA"), ids = c("a", "b", "c"))}
#' @importFrom seqinr translate
#' @importFrom muscle muscle
#' @export

# msa3 <- function(sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE){
#
#   tr <- function(seq, genetic_code = 1){
#     seq <- gsub(" ", "", seq)
#     seq <- strsplit(seq, "")[[1]]
#     output <- paste(translate(seq, numcode = genetic_code), collapse = "")
#     return(output)
#   }
#   if (length(sequences) < 2) {
#     stop("At least two sequences are required!")
#   } else if (length(sequences) != length(ids)) {
#     stop("The number of sequences and sequences' ids doesn't match!")
#   }
#
#   if (requireNamespace('Biostrings', quietly = TRUE)){
#     if (seqtype == "prot"){
#       seq <- Biostrings::AAStringSet(sequences)
#     } else if (seqtype == "dna"){
#       seq <- Biostrings::DNAStringSet(sequences)
#     }
#   } else {
#     stop("You must install the package Biostrings in order to use this function")
#   }
#   aln1 <- muscle(seq)
#   aln <- list()
#   aln$seq <- sequences
#   aln$ids <- ids
#   aln$aln <- as.character(aln1)
#   l <- sapply(aln$aln, function(x) strsplit(x, split = ""))
#   aln$ali <- matrix(unlist(l), nrow = length(sequences),
#                     byrow = TRUE)
#   if (sfile != FALSE) {
#     for (i in 1:length(aln$aln)) {
#       t <- paste(">", aln$ids[i], sep = "")
#       cat(t, file = sfile, append = TRUE)
#       if (seqtype == "cds"){
#         tt <- paste("\n", aln$cod[i], "\n", sep = "")
#       } else {
#         tt <- paste("\n", aln$aln[i], "\n", sep = "")
#       }
#       cat(tt, file = sfile, append = TRUE)
#     }
#   }
#   return(aln)
# }



