# =====================================================
# Global variables (avoid NOTE in R CMD check)
# =====================================================
utils::globalVariables(c(
  "data_df",
  "nt_pos1", "nt_pos2", "nt_pos3", "nt_pos4", "nt_pos5"
))

# =====================================================
# dna_encoding()
# =====================================================

#' @title Encode DNA Sequences into One-Hot Representation
#' @description Converts a vector of DNA strings into a one-hot encoded data frame.
#' Each position becomes a factor with levels A, T, C, G.
#'
#' @param dna_strings Character vector of DNA sequences (all same length).
#'
#' @return A data frame where each column corresponds to a nucleotide position,
#' and each cell contains a factor: A/T/C/G.
#'
#'
#' @export
#'
#' @examples
#' dna_encoding(c("ATCGA", "TGCAT"))
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn, byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

# =====================================================
# prediction_multiple()
# =====================================================

#' @title Predict m6A Modification Status for Multiple Sequences
#' @description Generates predictions for multiple RNA sequences using a trained
#' random forest model. Returns a data frame containing the input features plus:
#' - predicted_m6A_prob
#' - predicted_m6A_status
#'
#' @param ml_fit A trained randomForest model.
#' @param feature_df Data frame with required columns:
#' gc_content, RNA_type, RNA_region, exon_length, distance_to_junction,
#' evolutionary_conservation, DNA_5mer.
#' @param positive_threshold Numeric threshold for "Positive" label.
#'
#' @return Data frame with prediction results appended.
#'
#' @import randomForest
#' @export
#'
#' @examples
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' test_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                 package = "m6APrediction"))
#' head(prediction_multiple(rf_model, test_df))
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {

  # Required columns check
  required_cols <- c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  )

  stopifnot(all(required_cols %in% colnames(feature_df)))

  set.seed(123)

  # Ensure correct factor levels
  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = ml_fit$forest$xlevels$RNA_type
  )
  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = ml_fit$forest$xlevels$RNA_region
  )

  feature_df$DNA_5mer <- as.character(feature_df$DNA_5mer)

  # Add one-hot encoding of 5-mer
  dna_seq <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df, dna_seq)

  # Prediction
  preds <- predict(ml_fit, newdata = feature_df, type = "prob")
  predicted_m6A_prob <- preds[, "Positive"]

  predicted_m6A_status <- ifelse(
    predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative"
  )

  feature_df$predicted_m6A_prob <- predicted_m6A_prob
  feature_df$predicted_m6A_status <- predicted_m6A_status

  # Remove one-hot columns from output
  feature_df <- subset(
    feature_df,
    select = -c(nt_pos1, nt_pos2, nt_pos3, nt_pos4, nt_pos5)
  )

  return(feature_df)
}

# =====================================================
# prediction_single()
# =====================================================

#' @title Predict m6A Modification Status for a Single Sequence
#' @description Predicts m6A modification probability and status for a single RNA sequence.
#'
#' @param ml_fit Trained random forest model.
#' @param gc_content Numeric.
#' @param RNA_type Character or factor (will be converted).
#' @param RNA_region Character or factor (will be converted).
#' @param exon_length Numeric.
#' @param distance_to_junction Numeric.
#' @param evolutionary_conservation Numeric.
#' @param DNA_5mer Character 5-mer.
#' @param positive_threshold Numeric threshold.
#'
#' @return Named vector with:
#' - predicted_m6A_prob
#' - predicted_m6A_status
#'
#' @export
#'
#' @examples
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction"))
#'
#' prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATCGAT"
#' )
prediction_single <- function(
    ml_fit,
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer,
    positive_threshold = 0.5
) {
  set.seed(123)

  # Match factor levels
  RNA_type <- factor(RNA_type, levels = ml_fit$forest$xlevels$RNA_type)
  RNA_region <- factor(RNA_region, levels = ml_fit$forest$xlevels$RNA_region)

  # Construct single-row data frame
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  # Use prediction_multiple() to handle encoding and prediction
  results_df <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )

  out <- unlist(results_df[1, c("predicted_m6A_prob", "predicted_m6A_status")])
  return(out)
}
