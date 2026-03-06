SCRIPT_NAME = "MassSpectFPAutoTool for Mascot"
SCRIPT_VERSION = "1.0.0"
SCRIPT_AUTHOR = "Spike Murphy Müller"

# ---------------------------------
# USER INFO
# ---------------------------------


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
# MASCOT PARAMETERS
# ---------------------------------

mascot_parameters <- list(
  
  USERNAME = "USERNAME",
  USEREMAIL = "USEREMAIL",
  
  SEARCH = "PMF",
  FORMVER = 1.01,
  
  DB = "SwissProt", # database
  
  CLE = "Trypsin", # restriction enzyme
  PFA = 1, # max missed cleavage sites
  
  MODS = "", # definite modifications
  IT_MODS = "Oxidation (M)", # possible modifications
  
  TOL = 0.3, # tolerance
  TOLU = "Da", # toleracne unit
  
  MASS = "Monoisotopic", # Monoisotopic/Average Peaks
  CHARGE = "1+",
  
  REPORT = "AUTO" # amount of reports
)

# convert parameter list to readable text for log
param_lines <- c()
for (n in names(mascot_parameters)) {
  param_lines <- c(
    param_lines,
    paste(n, "=", mascot_parameters[[n]])
  )
}
param_lines <- c(param_lines, "")

# ---------------------------------
# SELECT FILES
# ---------------------------------

cat("\nSelect a spectrum file\n")

f <- file.choose()
files <- c(f)

cat("\n---------------------------------\n")
cat("FILE SELECTED\n")
cat("---------------------------------\n")
print(files)

# ---------------------------------
# ASK USER WHICH SEARCHES TO RUN
# ---------------------------------

ans_mascot <- readline("Run Mascot search? (y/n): ")
run_mascot <- tolower(ans_mascot) == "y"

cat("\nSearch settings:\n")
cat("Mascot:", run_mascot, "\n")

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
  "MassSpectFPAutoTool for Mascot Log",
  paste("Start time:", format(start_time,"%Y-%m-%d %H:%M:%S")),
  paste("Files processed:", length(files)),
  "=====================================",
  "",
  "-------------------------------------",
  "Mascot Parameters",
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
# MASCOT PMF SUBMISSION
# ---------------------------------

submit_mascot <- function(mono, mascot_parameters){
  
  if(nrow(mono) == 0) return(NULL)
  
  peaks <- mono$mz
  peak_text <- paste(peaks, collapse="\n")
  
  html <- paste0(
    '<!DOCTYPE html>
<html>
<body onload="document.forms[0].submit()">

<form method="POST"
action="https://www.matrixscience.com/cgi/nph-mascot.exe?1+-batch+-format+msr"
enctype="multipart/form-data">

<input type="hidden" name="USERNAME" value="', mascot_parameters$USERNAME, '">
<input type="hidden" name="USEREMAIL" value="', mascot_parameters$USEREMAIL, '">

<input type="hidden" name="SEARCH" value="', mascot_parameters$SEARCH, '">
<input type="hidden" name="FORMVER" value="', mascot_parameters$FORMVER, '">

<input type="hidden" name="DB" value="', mascot_parameters$DB, '">

<input type="hidden" name="CLE" value="', mascot_parameters$CLE, '">
<input type="hidden" name="PFA" value="', mascot_parameters$PFA, '">

<input type="hidden" name="MODS" value="', mascot_parameters$MODS, '">
<input type="hidden" name="IT_MODS" value="', mascot_parameters$IT_MODS, '">

<input type="hidden" name="TOL" value="', mascot_parameters$TOL, '">
<input type="hidden" name="TOLU" value="', mascot_parameters$TOLU, '">

<input type="hidden" name="MASS" value="', mascot_parameters$MASS, '">

<input type="hidden" name="CHARGE" value="', mascot_parameters$CHARGE, '">

<input type="hidden" name="REPORT" value="', mascot_parameters$REPORT, '">

<textarea name="QUE">', peak_text, '</textarea>

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
    if(run_mascot){
      submit_mascot(mono, mascot_parameters)
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

