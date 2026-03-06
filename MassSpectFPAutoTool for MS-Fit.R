SCRIPT_NAME = "MassSpectFPAutoTool for MS-Fit"
SCRIPT_VERSION = "1.0.0"
SCRIPT_AUTHOR = "Spike Murphy Müller"

# ---------------------------------
# USER INFO
# ---------------------------------

# in the MS-FIT FUNCTION update parameters if necessary
## (available parameters at: https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msfitupload)
## also run .random.concat library search if possible for FDR (false disccovery rate) calculation
## FDR = hits(library)/hits(library.random.concat)
# run this script with "command+shift+enter"
# select the first ms spectrum
# chose y/n for another spectrum or not
# when n start processing all chosen files one after the other
## 1. preprocess
## 2. try and select peaks and pic monoiotopic peaks with SNR 30, if not 20-40 monoiotopic peaks adjust SNR
## 3. submit ms-fit search
## 4. export peak list and png of peaks and monoisotopic peaks
# at the end export a log file

# ---------------------------------
# INSTALL LIBRARIES IF NEEDED
# ---------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("MALDIquant", quietly = TRUE))
  BiocManager::install("MALDIquant")

if (!requireNamespace("MALDIquantForeign", quietly = TRUE))
  BiocManager::install("MALDIquantForeign")

library(MALDIquant)
library(MALDIquantForeign)

# ---------------------------------
# MS-FIT PARAMETERS
# ---------------------------------

msfit_parameters <- list(
  
  database = "SwissProt.2025.02.26",
  dna_frame_translation = 3,
  n_term_aa_limit = "",
  species = "All",
  output_type = "HTML",
  sort_type = "Score Sort",
  ms_report_homologous_proteins = "Interesting",
  
  user_protein_sequence = "", # leave empty unless used
  
  # pre search parameters
  prot_low_mass = 1000,
  prot_high_mass = 125000,
  pi_low = 3.0,
  pi_high = 10.0,
  
  # digest parameters
  enzyme = "Trypsin",
  missed_cleavages = 1,
  
  # modifications
  const_mod = "",
  variable_mod = "Oxidation of M",
  
  # hit parameters
  parent_mass_convert = "monoisotopic",
  parent_mass_tolerance = 0.3,
  tolerance_units = "Da",
  
  min_matches = 4,
  min_parent_matches = 1,
  max_modifications = 1,
  max_reported_hits = 5,
  
  display_graph = 1,  # 0 = False, 1 = True
  mowse_on = 1, # 0 = False, 1 = True
  mowse_pfactor = 0.4,
  
  # instrument
  instrument = "MALDI-TOFTOF",
  data_format = "PP M/Z Charge"
  
)

# convert parameter list to readable text for log
param_lines <- c()
for (n in names(msfit_parameters)) {
  param_lines <- c(
    param_lines,
    paste(n, "=", msfit_parameters[[n]])
  )
}
param_lines <- c(param_lines, "")

# ---------------------------------
# SELECT FILES
# ---------------------------------

files <- c()

repeat {
  cat("\nSelect a spectrum file\n")
  f <- file.choose()
  cat("Selected file:\n")
  print(f)
  files <- c(files, f)
  ans <- readline("Add another file? (y/n): ")
  if (tolower(ans) != "y") break
}

cat("\n---------------------------------\n")
cat("FILES SELECTED\n")
cat("---------------------------------\n")
print(files)

# ---------------------------------
# ASK USER WHICH SEARCHES TO RUN
# ---------------------------------

ans_msfit <- readline("Run MS-Fit search? (y/n): ")
run_msfit <- tolower(ans_msfit) == "y"

ans_random <- readline("Run MS-Fit for FDR analysis? (y/n): ")
run_random <- tolower(ans_random) == "y"

cat("\nSearch settings:\n")
cat("MS-Fit:", run_msfit, "\n")
cat("MS-Fit .random.concat:", run_random, "\n\n")

# ---------------------------------
# CREATE RESULTS DIRECTORY
# ---------------------------------

base_dir <- dirname(files[1])
output_dir <- file.path(base_dir, "results")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ---------------------------------
# PREPARE LOG
# ---------------------------------

logfile <- file.path(output_dir, "processing_log.txt")
log_lines <- c()

start_time <- Sys.time()

log_lines <- c(
  log_lines,
  "=====================================",
  "MassSpectFPAutoTool for MS-Fit Log",
  paste("Start time:", format(start_time,"%Y-%m-%d %H:%M:%S")),
  paste("Files processed:", length(files)),
  "=====================================",
  "",
  "-------------------------------------",
  "MS-Fit Parameters",
  "-------------------------------------",
  param_lines,
  "-------------------------------------",
  "FILES",
  "-------------------------------------",
  paste(seq_along(files), basename(files)),
  ""
)

# summary storage
summary <- data.frame(
  File=character(),
  SNR=numeric(),
  Peaks=integer(),
  Duration=numeric(),
  stringsAsFactors=FALSE
)

# ---------------------------------
# MS-FIT FUNCTION
# ---------------------------------

submit_msfit <- function(mono, msfit_parameters, random_concat = FALSE){
  
  if(nrow(mono) == 0) return(NULL)
  
  peaks <- mono$mz
  peak_text <- paste(peaks, collapse="\n")
  
  db <- msfit_parameters$database
  if(random_concat){
    db <- paste0(db, ".random.concat")
  }
  
  html <- paste0(
    '<!DOCTYPE html>
<html>
<body onload="document.forms[0].submit()">

<form method="POST" action="https://prospector.ucsf.edu/prospector/cgi-bin/mssearch.cgi">

<input type="hidden" name="search_name" value="msfit">
<input type="hidden" name="report_title" value="MS-Fit">
<input type="hidden" name="version" value="6.8.1">
<input type="hidden" name="data_source" value="Data Paste Area">

<input type="hidden" name="ms_peak_exclusion" value="0">
<input type="hidden" name="ms_mass_exclusion" value="0">
<input type="hidden" name="ms_matrix_exclusion" value="0">
<input type="hidden" name="detailed_report" value="1">

<input type="hidden" name="database" value="', db, '">

<input type="hidden" name="species" value="', msfit_parameters$species, '">

<input type="hidden" name="dna_frame_translation" value="', msfit_parameters$dna_frame_translation, '">
<input type="hidden" name="n_term_aa_limit" value="', msfit_parameters$n_term_aa_limit, '">

<input type="hidden" name="output_type" value="', msfit_parameters$output_type, '">
<input type="hidden" name="sort_type" value="', msfit_parameters$sort_type, '">
<input type="hidden" name="ms_report_homologous_proteins" value="', msfit_parameters$ms_report_homologous_proteins, '">

<input type="hidden" name="enzyme" value="', msfit_parameters$enzyme, '">
<input type="hidden" name="missed_cleavages" value="', msfit_parameters$missed_cleavages, '">

<input type="hidden" name="const_mod" value="', msfit_parameters$const_mod, '">
<input type="hidden" name="mod_AA" value="', msfit_parameters$variable_mod, '">

<input type="hidden" name="parent_mass_convert" value="', msfit_parameters$parent_mass_convert, '">
<input type="hidden" name="ms_parent_mass_tolerance" value="', msfit_parameters$parent_mass_tolerance, '">
<input type="hidden" name="ms_parent_mass_tolerance_units" value="', msfit_parameters$tolerance_units, '">
<input type="hidden" name="ms_parent_mass_systematic_error" value="0">

<input type="hidden" name="min_matches" value="', msfit_parameters$min_matches, '">
<input type="hidden" name="min_parent_ion_matches" value="', msfit_parameters$min_parent_matches, '">
<input type="hidden" name="ms_max_modifications" value="', msfit_parameters$max_modifications, '">
<input type="hidden" name="ms_max_reported_hits" value="', msfit_parameters$max_reported_hits, '">

<input type="hidden" name="mowse_on" value="', msfit_parameters$mowse_on, '">
<input type="hidden" name="mowse_pfactor" value="', msfit_parameters$mowse_pfactor, '">

<input type="hidden" name="ms_prot_low_mass" value="', msfit_parameters$prot_low_mass, '">
<input type="hidden" name="ms_prot_high_mass" value="', msfit_parameters$prot_high_mass, '">

<input type="hidden" name="low_pi" value="', msfit_parameters$pi_low, '">
<input type="hidden" name="high_pi" value="', msfit_parameters$pi_high, '">
<input type="hidden" name="full_pi_range" value="1">

<input type="hidden" name="instrument_name" value="', msfit_parameters$instrument, '">
<input type="hidden" name="data_format" value="', msfit_parameters$data_format, '">

<input type="hidden" name="display_graph" value="', msfit_parameters$display_graph, '">
<input type="hidden" name="comment" value="">

<input type="hidden" name="output_filename" value="results">

<textarea name="user_protein_sequence" rows="6" cols="70">', msfit_parameters$user_protein_sequence, '</textarea>
<textarea name="data">', peak_text, '</textarea>

</form>
</body>
</html>'
  )
  
  tmp <- tempfile(fileext = ".html")
  writeLines(html, tmp)
  
  browseURL(tmp)
}

# ---------------------------------
# PROCESS FILES
# ---------------------------------

for (i in seq_along(files)) {
  f <- files[i]
  file_start <- Sys.time()
  
  cat("\nProcessing:", basename(f), "\n")
  
  try({
    # preprocessing of spectrums
    spectra_raw <- import(f)
    
    spectra_processed <- spectra_raw
    spectra_processed <- transformIntensity(spectra_processed, method="sqrt")
    spectra_processed <- smoothIntensity(spectra_processed, method="SavitzkyGolay", halfWindowSize=10)
    spectra_processed <- removeBaseline(spectra_processed, method="SNIP", iterations=100)
    spectra_processed <- calibrateIntensity(spectra_processed, method="TIC")
    
    spectra <- spectra_processed
    
    # initilize variables
    snr <- 30
    attempts <- 0
    
    repeat {
      
      peaks <- detectPeaks(
        spectra,
        SNR = snr,
        halfWindowSize = 20
      )
      
      # remove trypsin peaks
      trypsin_peaks <- c(
        842.51, 
        1045.56, 
        2211.10
      )
      tol <- 3
      mz_vals <- mass(peaks[[1]])
      remove_idx <- sapply(mz_vals, function(m)
        any(abs(m - trypsin_peaks) < tol)
      )
      trypsin_removed <- data.frame(
        mz = mz_vals[remove_idx],
        intensity = intensity(peaks[[1]])[remove_idx]
      )
      # create monoisotopic trypsin peaks
      trypsin_mono <- trypsin_removed[order(trypsin_removed$mz), ]
      if(nrow(trypsin_mono) > 0){
        trypsin_mono <- trypsin_mono[c(TRUE, diff(trypsin_mono$mz) > 1.2), ]
      }
      peaks[[1]] <- peaks[[1]][!remove_idx]
      
      # remove keratin peaks
      keratin_peaks <- c(
        1165.59,
        1191.64,
        1433.74,
        1479.76,
        1783.88,
        1905.95,
        1979.96,
        2045.10,
        2155.10
      )
      tol <- 3
      mz_vals <- mass(peaks[[1]])
      keratin_idx <- sapply(mz_vals, function(m)
        any(abs(m - keratin_peaks) < tol)
      )
      keratin_removed <- data.frame(
        mz = mz_vals[keratin_idx],
        intensity = intensity(peaks[[1]])[keratin_idx]
      )
      # create monoisotopic keratin peaks
      keratin_mono <- keratin_removed[order(keratin_removed$mz), ]
      if(nrow(keratin_mono) > 0){
        keratin_mono <- keratin_mono[c(TRUE, diff(keratin_mono$mz) > 1.2), ]
      }
      peaks[[1]] <- peaks[[1]][!keratin_idx]
      
      p <- peaks[[1]]
      peaklist <- data.frame(
        mz = mass(p),
        intensity = intensity(p)
      )
      
      # remove isotopic peaks 
      mono <- peaklist[order(peaklist$mz), ]
      mono <- mono[c(TRUE, diff(mono$mz) > 1.2), ]
      
      # check amout of peaks and adjust SNR
      n_peaks <- nrow(mono)
      if (n_peaks >= 20 && n_peaks <= 40) break
      attempts <- attempts + 1
      if (n_peaks < 20) {
        snr <- snr * 0.8
      } else {
        snr <- snr * 1.2
      }
      if (attempts > 10) break
    }
    
    # ---------------------------------
    # RUN SEARCH
    # ---------------------------------
    
    # submit searches depending on user choice
    if(nrow(mono) > 0){
      if(run_msfit){
        submit_msfit(mono, msfit_parameters)
      }
      if(run_random){
        submit_msfit(mono, msfit_parameters, random_concat = TRUE)
      }
    }
    
    # ---------------------------------
    # Base Filename
    # ---------------------------------
    
    base <- tools::file_path_sans_ext(basename(f))
    
    # ---------------------------------
    # EXPORT original spectrum
    # ---------------------------------
    
    png(
      file.path(output_dir,
                paste0(base,"_spectrum_original.png")),
      width = 1200,
      height = 800
    )
    
    plot(
      spectra_raw[[1]],
      main = paste(
        "Original Spectrum |",
        basename(f)
      ),
      sub = ""
    )
    
    dev.off()
    
    # ---------------------------------
    # EXPORT preprocessed spectrum
    # ---------------------------------
    
    png(
      file.path(output_dir,
                paste0(base,"_spectrum_preprocessed.png")),
      width = 1200,
      height = 800
    )
    
    plot(
      spectra_processed[[1]],
      main = paste(
        "Preprocessed Spectrum |",
        basename(f)
      ),
      sub = ""
    )
    
    dev.off()
    
    # ---------------------------------
    # EXPORT DETECTED PEAKS PLOT
    # ---------------------------------
    
    png(
      file.path(output_dir,
                paste0(base,"_spectrum_peaks.png")),
      width=1200,
      height=800
    )
    
    plot(
      spectra[[1]],
      main=paste(
        "Detected Peaks |",
        basename(f),
        "| SNR =",round(snr,2)
      ),
      sub=""
    )
    
    points(peaks[[1]],col="blue",pch=19)
    
    points(
      trypsin_removed$mz,
      trypsin_removed$intensity,
      col="orange",
      pch=19,
    )
    
    points(
      keratin_removed$mz,
      keratin_removed$intensity,
      col="brown",
      pch=19,
    )
    
    dev.off()
    
    # ---------------------------------
    # EXPORT MONOISOTOPIC PEAKS PLOT
    # ---------------------------------
    
    png(
      file.path(output_dir,
                paste0(base,"_spectrum_monoisotopic_peaks.png")),
      width=1200,
      height=800
    )
    
    plot(
      spectra[[1]],
      main=paste(
        "Monoisotopic Peaks |",
        basename(f)
      ),
      sub=""
    )
    
    points(
      mono$mz,
      mono$intensity,
      col="forestgreen",
      pch=19
    )
    
    points(
      trypsin_mono$mz,
      trypsin_mono$intensity,
      col="orange",
      pch=19,
    )
    
    points(
      keratin_mono$mz,
      keratin_mono$intensity,
      col="brown",
      pch=19,
    )
    
    dev.off()
    
    # ---------------------------------
    # EXPORT PEAK LIST
    # ---------------------------------
    
    outfile <- file.path(
      output_dir,
      paste0(base,"_monoisotopic_peaks.txt")
    )
    
    write.table(
      mono$mz,
      outfile,
      row.names=FALSE,
      col.names=FALSE
    )
    
    # ---------------------------------
    # LOG
    # ---------------------------------
    
    file_end <- Sys.time()
    
    duration <- round(
      as.numeric(difftime(file_end,file_start,units="secs")),2
    )
    
    # add to summary
    summary <- rbind(
      summary,
      data.frame(
        File=base,
        SNR=round(snr,2),
        Peaks=n_peaks,
        Duration=duration
      )
    )
    
    # add log entry
    log_lines <- c(
      log_lines,
      paste0("Processing file ",i,"/",length(files),": ",basename(f)),
      paste("Start:",format(file_start,"%H:%M:%S")),
      paste("Mono peaks detected:",n_peaks),
      paste("SNR adjustments:",attempts),
      paste("Final SNR:",round(snr,2)),
      paste("Duration:",duration,"sec"),
      paste("Output:",basename(outfile)),
      ""
    )
    
  }, silent=TRUE)
}

# ---------------------------------
# SUMMARY
# ---------------------------------

log_lines <- c(
  log_lines,
  "-------------------------------------",
  "SUMMARY",
  "-------------------------------------",
  sprintf("%-8s %-8s %-8s %-12s",
          "File","SNR","Peaks","Duration")
)

for(i in 1:nrow(summary)){
  
  log_lines <- c(
    log_lines,
    sprintf("%-8s %-8s %-8s %-12s",
            summary$File[i],
            summary$SNR[i],
            summary$Peaks[i],
            summary$Duration[i])
  )
}

end_time <- Sys.time()

runtime <- round(
  as.numeric(difftime(end_time,start_time,units="secs")),2
)

log_lines <- c(
  log_lines,
  "",
  paste("Total runtime:",runtime,"sec"),
  paste("End time:",format(end_time,"%Y-%m-%d %H:%M:%S")),
  "====================================="
)

# ---------------------------------
# WRITE LOG FILE
# ---------------------------------

writeLines(log_lines, logfile)

cat("\nProcessing finished.\n")
cat("Log saved to:\n")
print(logfile)

