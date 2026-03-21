# =================================================== #
# ===== MAIN FUNCTION =============================== #
# =================================================== #

run_processing <- function (FILE, CLEAVAGE = "Trypsin", KERATIN = TRUE, TAG = NULL, EXCLUDE = c(), SNR = 30, PEAKS = c(30, 10), PLOTS = TRUE, TRANSFORM = "sqrt", SMOOTH = c("SavitzkyGolay", 10), BASELINE = c("SNIP", 100), CALIBRATE = "TIC"){
  # read file(s)
  raw_spectrum <- MALDIquantForeign::import(FILE)
  
  # preprossessing
  preprocessed_spectrum <- raw_spectrum |>
    MALDIquant::transformIntensity(method = TRANSFORM) |> # Reduce dominance of very large peaks with square root scaling
    MALDIquant::smoothIntensity(method = SMOOTH[1], halfWindowSize = as.numeric(SMOOTH[2])) |> # Reduce noise while preserving peak shape
    MALDIquant::removeBaseline(method = BASELINE[1], iterations = as.numeric(BASELINE[2])) |> # Remove baseline with iterative clipping algorithm
    MALDIquant::calibrateIntensity(method = CALIBRATE) # Calibrate with total ion current normalization
  
  # processing
  ## initialize
  sigal_noise_ratio <- SNR
  number_peaks <- 0
  min_peaks = PEAKS[1] - PEAKS[2]
  max_peaks = PEAKS[1] + PEAKS[2]
  
  # source: https://home.pavlab.msl.ubc.ca/wp-content/uploads/2013/06/Notes-on-trouble-shooting-LCMS-contamination-full.pdf
  if (CLEAVAGE == "Trypsin"){
    trypsin_peaks <- c(802.4, 842.87, 1046.00, 1567.6, 1713.8, 1940.9, 2083.4, 
                       2211.10)
  }

  
  # source: https://home.pavlab.msl.ubc.ca/wp-content/uploads/2013/06/Notes-on-trouble-shooting-LCMS-contamination-full.pdf
  if (KERATIN == TRUE){
    keratin_peaks <- c(704.4, 809.4, 827.4, 875.0, 973.5, 995.5, 1000.6, 1003.5, 
                       1031.6, 1033.5, 1036.5, 1060.6, 1066.0, 1129.6, 1141.5, 
                       1165.6, 1179.6, 1190.6, 1262.6, 1265.6, 1277.7, 1278.5, 
                       1301.7, 1302.7, 1307.7, 1315.7, 1340.6, 1344.7, 1350.7, 
                       1357.7, 1371.7, 1381.6, 1383.7, 1390.7, 1394.56, 1418.7, 
                       1419.7, 1442.8, 1453.8, 1475.8, 1476.7, 1493.7, 1549.6, 
                       1586.8, 1599.8, 1657.8, 1708.8, 1716.8, 1765.7, 1791.7, 
                       1792.9, 1798.1, 1837.9, 1847.8, 1851.9, 1994.0, 2082.9, 
                       2109.0, 2171.0, 2184.1, 2240.1, 2330.49, 2501.2, 2510.1, 
                       2565.88, 2581.1, 2705.1, 2746.4, 2872.4, 2902.7, 2904.4, 
                       2932.5, 3223.2)
  }
  if (!is.null(TAG)) {
    if (TAG == "GFP") {
      # source: sequence from uniprot, peaks from MS-Digest (https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msdigest)
      tag_peaks <- c(806.3502, 821.3941, 919.536, 968.4432, 982.4952, 
                     1050.5214, 1062.6194, 1224.7021, 1266.5783, 1282.5732, 
                     1347.6579, 1477.7645, 1503.6598, 1533.8271, 1542.7911, 
                     1592.7315, 1608.7264, 1902.9418, 1918.9368, 1958.9706, 
                     1973.9062, 1989.9011, 2068.9545, 2084.9495, 2230.0597, 
                     2246.0546, 2398.2264, 2437.2609, 2590.268, 2606.2629, 
                     2622.2578, 2746.3691, 2762.364, 2778.3589, 2783.4284, 
                     2799.4233, 2937.3836, 2953.3785, 3148.5998, 3169.5638, 
                     3185.5587, 3921.9029)
    } else if (TAG == "RFP"){
      # source: sequence from uniprot, peaks from MS-Digest (https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msdigest)
      tag_peaks <- c(857.4338, 873.4287, 895.536, 924.4825, 1036.5608, 
                     1052.5557, 1052.5775, 1055.5156, 1060.615, 1095.6157, 
                     1109.5851, 1164.5684, 1167.5099, 1183.5048, 1183.6106, 
                     1295.6049, 1311.5998, 1395.7056, 1405.7474, 1544.8108, 
                     1635.8199, 1651.8149, 1656.8632, 1696.8112, 1712.8061, 
                     1763.9149, 1779.9098, 1857.9204, 1862.9833, 1873.9153, 
                     1878.9782, 2208.1124, 2225.0357, 2350.2152, 2366.2101, 
                     2536.3235, 2668.2485, 2697.4075, 2852.3117, 2868.3066, 
                     3079.4751, 3095.47, 3102.3858, 3118.3807, 3323.5347, 
                     3339.5297, 3343.5648, 3357.5553, 3359.5597, 3373.5502, 
                     3862.881)
    }
  }
  
  if (length(EXCLUDE) > 0){
    custom_peaks <- EXCLUDE
  }
  
  ## loop
  while (number_peaks < min_peaks || number_peaks > max_peaks){
    ### detect peaks
    peaks <- MALDIquant::detectPeaks(
      preprocessed_spectrum,
      SNR = sigal_noise_ratio,
      halfWindowSize = 20 # ?????
    )
    ### remove contaminants (protease, keratin, tag)
    peaks_mz <- MALDIquant::mass(peaks[[1]])
    peaks_int <- MALDIquant::intensity(peaks[[1]])
    #### detect matches
    tolerance <- 2
    
    cleavage_counter <- 0
    keratin_counter <- 0
    tag_counter <- 0
    custom_counter <- 0
    
    keep <- rep(TRUE, length(peaks_mz))
    keep_cleavage <- rep(FALSE, length(peaks_mz))
    keep_keratin <- rep(FALSE, length(peaks_mz))
    keep_tag <- rep(FALSE, length(peaks_mz))
    keep_custom <- rep(FALSE, length(peaks_mz))
    
    if (CLEAVAGE == "Trypsin"){
      trypsin_matches <- c()
      for (mz in peaks_mz){
        cleavage_counter <- cleavage_counter + 1
        for (peak in trypsin_peaks) {
          distance <- abs(mz - peak)
          if (distance < tolerance) {
            trypsin_matches <- c(trypsin_matches, mz)
            keep[cleavage_counter] <- FALSE
            keep_cleavage[cleavage_counter] <- TRUE
            break
          }
        }
      }
    }
    
    if (KERATIN == TRUE){
      keratin_matches <- c()
      for (mz in peaks_mz){
        keratin_counter <- keratin_counter + 1
        for (peak in keratin_peaks) {
          distance <- abs(mz - peak)
          if (distance < tolerance) {
            keratin_matches <- c(keratin_matches, mz)
            keep[keratin_counter] <- FALSE
            keep_keratin[keratin_counter] <- TRUE
            break
          }
        }
      }
    }
    
    if (!is.null(TAG)){
      tag_matches <- c()
      for (mz in peaks_mz){
        tag_counter <- tag_counter + 1
        for (peak in tag_peaks) {
          distance <- abs(mz - peak)
          if (distance < tolerance) {
            tag_matches <- c(tag_matches, mz)
            keep[tag_counter] <- FALSE
            keep_tag[tag_counter] <- TRUE
            break
          }
        }
      }
    }
    
    if (length(EXCLUDE) > 0){
      custom_matches <- c()
      for (mz in peaks_mz){
        custom_counter <- custom_counter + 1
        for (peak in custom_peaks) {
          distance <- abs(mz - peak)
          if (distance < tolerance) {
            custom_matches <- c(custom_matches, mz)
            keep[custom_counter] <- FALSE
            keep_custom[custom_counter] <- TRUE
            break
          }
        }
      }
    }
    
    ### remove contaminant peaks
    peaks_mz_clean <- peaks_mz[keep]
    peaks_mz_clean_cleavage <- peaks_mz[keep_cleavage]
    peaks_mz_clean_keratin <- peaks_mz[keep_keratin]
    peaks_mz_clean_tag <- peaks_mz[keep_tag]
    peaks_mz_clean_custom <- peaks_mz[keep_custom]
    peaks_int_clean <- peaks_int[keep]
    peaks_int_clean_cleavage <- peaks_int[keep_cleavage]
    peaks_int_clean_keratin <- peaks_int[keep_keratin]
    peaks_int_clean_tag <- peaks_int[keep_tag]
    peaks_int_clean_custom <- peaks_int[keep_custom]
    
    ### select monoisotopic peaks
    #### initiate
    isotope_distance <- c(0.8, 1.2)

    mono_counter <- 0
    isotope_counter <- 0
    
    keep_mono <- rep(TRUE, length(peaks_mz_clean))
    
    #### start loop
    for (mz1 in peaks_mz_clean){
      #### increment counter for proposed monoisotopic peak
      mono_counter <- mono_counter + 1
      #### initialize counter for proposed isotopic peak
      isotope_counter <- 0
      #### set the prior isotope first to the current proposed monoisotopic peak
      prior_isotope <- mz1
      for (mz2 in peaks_mz_clean){
        #### increment counter for proposed isotopic peak 
        isotope_counter <- isotope_counter + 1
        #### only continue once proposed isotopic peaks greater than proposed monoisotopic peak are reached
        if (mz2 > mz1){
          #### calculate the difference between current proposed isotopic peak value and current proposed monoisotopic peak value or last isotopic peak value
          difference = mz2 - prior_isotope
          #### continue if proposed isotopic peak was detected as isotope
          if (difference >= isotope_distance[1] && difference <= isotope_distance[2]){
            #### update logical vector
            keep_mono[isotope_counter] <- FALSE
            #### update  the prior isotope to the current detected isotopic peak
            prior_isotope <- mz2
          } else if (difference > isotope_distance[2]){
            #### break if no isotopic peak was detected within the specified tolerance
            break
          }
        }
      }
    }
    
    ### remove isotopic peaks
    peaks_mz_clean_mono <- peaks_mz_clean[keep_mono]
    peaks_int_clean_mono <- peaks_int_clean[keep_mono]
    
    ### check number of peak 
    number_peaks = length(peaks_mz_clean_mono)
    
    ### adjust signal_noise_ratio
    if (number_peaks < min_peaks){
      peak_differnce <- min_peaks - number_peaks
      sigal_noise_ratio <- sigal_noise_ratio - 1 * sqrt(peak_differnce)
    }else if (number_peaks > max_peaks){
      peak_differnce <- number_peaks - max_peaks
      sigal_noise_ratio <- sigal_noise_ratio + 1 * sqrt(peak_differnce)
    }
    
    ### generate list of monoisotopic peaks for contaminants (for comments see above)
    mono_counter <- 0
    
    keep_cleavage_mono <- rep(TRUE, length(peaks_mz_clean_cleavage))
    keep_keratin_mono <- rep(TRUE, length(peaks_mz_clean_keratin))
    keep_tag_mono <- rep(TRUE, length(peaks_mz_clean_tag))
    keep_custom_mono <- rep(TRUE, length(peaks_mz_clean_custom))
    
    #### cleavage
    if (number_peaks >= min_peaks && number_peaks <= max_peaks){
      for (mz1 in peaks_mz_clean_cleavage){
        mono_counter <- mono_counter + 1
        isotope_counter <- 0
        prior_isotope <- mz1
        for (mz2 in peaks_mz_clean_cleavage){
          isotope_counter <- isotope_counter + 1
          if (mz2 > mz1){
            difference = mz2 - prior_isotope
            if (difference >= isotope_distance[1] && difference <= isotope_distance[2]){
              keep_cleavage_mono[isotope_counter] <- FALSE
              prior_isotope <- mz2
            } else if (difference > isotope_distance[2]){
              break
            }
          }
        }
      }
    }
    #### keratin
    if (number_peaks >= min_peaks && number_peaks <= max_peaks){
      for (mz1 in peaks_mz_clean_keratin){
        mono_counter <- mono_counter + 1
        isotope_counter <- 0
        prior_isotope <- mz1
        for (mz2 in peaks_mz_clean_keratin){
          isotope_counter <- isotope_counter + 1
          if (mz2 > mz1){
            difference = mz2 - prior_isotope
            if (difference >= isotope_distance[1] && difference <= isotope_distance[2]){
              keep_keratin_mono[isotope_counter] <- FALSE
              prior_isotope <- mz2
            } else if (difference > isotope_distance[2]){
              break
            }
          }
        }
      }
    }
    #### tag
    if (number_peaks >= min_peaks && number_peaks <= max_peaks){
      for (mz1 in peaks_mz_clean_tag){
        mono_counter <- mono_counter + 1
        isotope_counter <- 0
        prior_isotope <- mz1
        for (mz2 in peaks_mz_clean_tag){
          isotope_counter <- isotope_counter + 1
          if (mz2 > mz1){
            difference = mz2 - prior_isotope
            if (difference >= isotope_distance[1] && difference <= isotope_distance[2]){
              keep_tag_mono[isotope_counter] <- FALSE
              prior_isotope <- mz2
            } else if (difference > isotope_distance[2]){
              break
            }
          }
        }
      }
    }
    #### custom
    if (number_peaks >= min_peaks && number_peaks <= max_peaks){
      for (mz1 in peaks_mz_clean_custom){
        mono_counter <- mono_counter + 1
        isotope_counter <- 0
        prior_isotope <- mz1
        for (mz2 in peaks_mz_clean_custom){
          isotope_counter <- isotope_counter + 1
          if (mz2 > mz1){
            difference = mz2 - prior_isotope
            if (difference >= isotope_distance[1] && difference <= isotope_distance[2]){
              keep_custom_mono[isotope_counter] <- FALSE
              prior_isotope <- mz2
            } else if (difference > isotope_distance[2]){
              break
            }
          }
        }
      }
    }
    #### remove isotopic peaks
    peaks_mz_clean_cleavage_mono <- peaks_mz_clean_cleavage[keep_cleavage_mono]
    peaks_int_clean_cleavage_mono <- peaks_int_clean_cleavage[keep_cleavage_mono]
    peaks_mz_clean_keratin_mono <- peaks_mz_clean_keratin[keep_keratin_mono]
    peaks_int_clean_keratin_mono <- peaks_int_clean_keratin[keep_keratin_mono]
    peaks_mz_clean_tag_mono <- peaks_mz_clean_tag[keep_tag_mono]
    peaks_int_clean_tag_mono <- peaks_int_clean_tag[keep_tag_mono]
    peaks_mz_clean_custom_mono <- peaks_mz_clean_custom[keep_custom_mono]
    peaks_int_clean_custom_mono <- peaks_int_clean_custom[keep_custom_mono]
    
  }
  
  # Plots
  ## Plot Raw Data
  plot_raw_data <- data.frame(mz = MALDIquant::mass(raw_spectrum[[1]]), intensity = MALDIquant::intensity(raw_spectrum[[1]]))
  
  plot_raw <- ggplot2::ggplot(data = plot_raw_data, ggplot2::aes(x = mz, y =  intensity)) + 
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(limits = c(0, max(MALDIquant::intensity(raw_spectrum[[1]])))) +
    ggplot2::labs(
      x = "m/z",
      y = "intensity",
      title = "Raw Spectrum"
    ) +
    ggplot2::theme_classic()
  
  ## Plot Preprocessed Data
  plot_preprocessed_data <- data.frame(mz = MALDIquant::mass(preprocessed_spectrum[[1]]), intensity = MALDIquant::intensity(preprocessed_spectrum[[1]]))

  plot_preprocessed <- ggplot2::ggplot(data = plot_preprocessed_data, ggplot2::aes(x = mz, y =  intensity)) + 
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(limits = c(0, max(MALDIquant::intensity(preprocessed_spectrum[[1]])))) +
    ggplot2::labs(
      x = "m/z",
      y = "intensity",
      title = "Preprocessed Spectrum"
    ) +
    ggplot2::theme_classic()
  
  ## Plot Detected Peaks
  plot_peaks_data <- data.frame(peaks = peaks_mz_clean, intensity = peaks_int_clean)
  
  plot_peaks <- ggplot2::ggplot(data = plot_preprocessed_data, ggplot2::aes(x = mz, y =  intensity)) + 
    ggplot2::geom_line() +
    ggplot2::geom_point(data = plot_peaks_data, ggplot2::aes(x = peaks, y = intensity), color = "lightblue4") +
    ggplot2::scale_y_continuous(limits = c(0, max(MALDIquant::intensity(preprocessed_spectrum[[1]])))) +
    ggplot2::labs(
      x = "m/z",
      y = "intensity",
      title = "Preprocessed Spectrum"
    ) +
    ggplot2::theme_classic()
  
  ## Plot Detected Peaks Clean Mono
  plot_peaks_mono_data <- data.frame(peaks = peaks_mz_clean_mono, intensity = peaks_int_clean_mono)
  
  plot_peaks_mono <- ggplot2::ggplot(data = plot_preprocessed_data, ggplot2::aes(x = mz, y =  intensity)) + 
    ggplot2::geom_line() +
    ggplot2::geom_point(data = plot_peaks_mono_data, ggplot2::aes(x = peaks, y = intensity), color = "red") +
    ggplot2::scale_y_continuous(limits = c(0, max(MALDIquant::intensity(preprocessed_spectrum[[1]])))) +
    ggplot2::labs(
      x = "m/z",
      y = "intensity",
      title = "Detected Monoisotopic Peaks For the Target Protein "
    ) +
    ggplot2::theme_classic()
  
  ## Plot Detected Peaks Clean Mono and Contaminants 
  plot_peaks_cleavage_mono <- data.frame(peaks = peaks_mz_clean_cleavage_mono, intensity = peaks_int_clean_cleavage_mono)
  plot_peaks_keratin_mono <- data.frame(peaks = peaks_mz_clean_keratin_mono, intensity = peaks_int_clean_keratin_mono)
  plot_peaks_tag_mono <- data.frame(peaks = peaks_mz_clean_tag_mono, intensity = peaks_int_clean_tag_mono)
  plot_peaks_custom_mono <- data.frame(peaks = peaks_mz_clean_custom_mono, intensity = peaks_int_clean_custom_mono)
  
  plot_peaks_mono_contaminants <- ggplot2::ggplot(data = plot_preprocessed_data, ggplot2::aes(x = mz, y =  intensity)) + 
    ggplot2::geom_line() +
    ggplot2::geom_point(data = plot_peaks_mono_data, ggplot2::aes(x = peaks, y = intensity, color = "Target")) +
    ggplot2::geom_point(data = plot_peaks_cleavage_mono, ggplot2::aes(x = peaks, y = intensity, color = "Cleavage")) +
    ggplot2::geom_point(data = plot_peaks_keratin_mono, ggplot2::aes(x = peaks, y = intensity, color = "Keratin")) +
    ggplot2::geom_point(data = plot_peaks_tag_mono, ggplot2::aes(x = peaks, y = intensity, color = "Tag")) +
    ggplot2::geom_point(data = plot_peaks_custom_mono, ggplot2::aes(x = peaks, y = intensity, color = "Custom")) +
    ggplot2::scale_color_manual(
      name = "Peak Type",
      values = c(
        "Target" = "red",
        "Cleavage" = "yellow",
        "Keratin" = "pink",
        "Tag" = "forestgreen",
        "Custom" = "purple"
      )
    ) +
    ggplot2::scale_y_continuous(limits = c(0, max(MALDIquant::intensity(preprocessed_spectrum[[1]])))) +
    ggplot2::labs(
      x = "m/z",
      y = "intensity",
      title = "Detected Monoisotopic Peaks For the Target Protein and the Removed Contaminants"
    ) +
    ggplot2::theme_classic()
  
  # return
  # TODO: change returnvalue according to agruments
  if (PLOTS == TRUE) {
    return(
      list(
        monoisotopic_peaks = plot_peaks_mono_data,
        monoisotopic_peaks_cleavage = plot_peaks_cleavage_mono,
        monoisotopic_peaks_keratin = plot_peaks_keratin_mono,
        monoisotopic_peaks_tag = plot_peaks_tag_mono,
        monoisotopic_peaks_custom = plot_peaks_custom_mono,
        plot_raw = plot_raw,
        plot_preprocessed = plot_preprocessed,
        plot_peaks = plot_peaks,
        plot_peaks_mono = plot_peaks_mono,
        plot_peaks_mono_contaminants = plot_peaks_mono_contaminants
      )
    )
  } else {
      return(
        list(
          monoisotopic_peaks = plot_peaks_mono_data,
          monoisotopic_peaks_cleavage = plot_peaks_cleavage_mono,
          monoisotopic_peaks_keratin = plot_peaks_keratin_mono,
          monoisotopic_peaks_tag = plot_peaks_tag_mono,
          monoisotopic_peaks_custom = plot_peaks_custom_mono
        )
      )
  }
  
}


# =================================================== #
# ===== ASSISTING FUNCTIONS ========================= #
# =================================================== #

processing_check_file <- function(FILE){
  if (length(FILE) != 1){
    stop("Argument `FILE` should consist of one file path.")
  }
    
  # check path
  if (!file.exists(FILE)){
    stop("Argument `FILE` indicates a non existing file.")
  }
  
  # supported extensions
  supported <- c("mzXML", "mzML", "imzML")
  
  # get extenstion
  split <- strsplit(FILE, "\\.")
  extension <- split[[1]][length(split[[1]])]
  
  # check extension
  if (!(extension %in% supported)){
    stop("Argument `FILE` has an unsupported file extension.")
  }
}


processing_check_cleavage <- function(CLEAVAGE){
  enzymes <- c("Trypsin")
  
  if (length(CLEAVAGE) != 1){
    warning("Argument `CLEAVAGE` should consist of one restriction enzyme.")
  } else if(!(CLEAVAGE %in% enzymes)){
    warning("Argument `CLEAVAGE` contains unexpected enzyme.")
  }
}


processing_check_contaminants <- function(KERATIN, TAG, EXCLUDE){
  # KERATIN
  if (length(KERATIN) != 1){
    warning("Argument `KERATIN` should contain one logical argument.")
  } else if (is.na(KERATIN)){
    warning("Argument `KERATIN` should be 'TRUE' or 'FALSE'.")
  }else if (KERATIN != TRUE && KERATIN != FALSE){
    warning("Argument `KERATIN` should be 'TRUE' or 'FALSE'.")
  }
    
  # TAG
  tags <- c("GFP", "RFP")
  
  if (!is.null(TAG)){
    if (length(TAG) != 1){
      warning("Argument `TAG` should consist of one restriction enzyme.")
    } else if(!(TAG %in% tags)){
      warning("Argument `TAG` contains unexpected enzyme. Try to generate a peak list for your tag with `search_msdigest()` and include the peaklist in the argument `EXCLUDE`")
    }
  }
  
  # EXCLUDE
  if (length(EXCLUDE) != 0){
    if (!is.numeric(EXCLUDE)){
      warning("Argument `EXCLUDE` should be a list of m/z values (peaks).")
    }
  }
  
}


processing_check_prerequisites <- function(SNR, PEAKS){
  
}


processing_check_plots <- function(PLOTS){
  
}


processing_check_processing <- function(TRANSFORM, SMOOTH, BASELINE, CALIBRATE){
  
}
