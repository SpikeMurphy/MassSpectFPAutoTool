# TODO: go through website HTML and orter arguments and local html

# =================================================== #
# ===== MAIN FUNCTION =============================== #
# =================================================== #

search_msfit <- function(
    PEAKS,
    SEQUENCE = NULL,
    
    DATABASE = "SwissProt.2025.02.26",
    SPECIES = "All",
    
    FRAME = 3,
    AALIMIT = NULL,

    OUTPUT = "HTML",
    SORT = "Score Sort",
    HOMOLOGUES = "Interesting",
    
    LOWMASS = 1000,
    HIGHMASS = 125000,
    LOWPI = 3.0,
    HIGHPI = 10.0,
    RANGEPI = "1",
    
    # digest parameters
    CLEAVAGE = "Trypsin",
    MISSEDC = 1,
    
    # modifications
    CONSTMODS = NULL,
    VARMODS = "Oxidation of M",
    
    # hit parameters
    MASS = "monoisotopic",
    TOL = 0.3,
    TOLU = "Da",
    
    MINMATCH = 4,
    MINPMATCH = 1,
    MAXMODS = 1,
    MAXHITS = 5,
    
    GRAPH = 1,  # 0 = False, 1 = True
    MOWSE = 1, # 0 = False, 1 = True
    MOWSEP = 0.4,
    
    INST = "MALDI-TOFTOF",
    FORMAT = "PP M/Z Charge",
    
    TEST = FALSE
  ){
  
  #check inputs
  msfit_check_peaks(PEAKS)
  msfit_check_sequence(SEQUENCE)
  msfit_check_database(DATABASE, SPECIES)
  msfit_check_frame(FRAME, AALIMIT)
  msfit_check_output(OUTPUT, SORT, HOMOLOGUES)
  msfit_check_prerequisites(LOWMASS, HIGHMASS, LOWPI, HIGHPI, RANGEPI)
  msfit_check_digest(CLEAVAGE, MISSEDC)
  msfit_check_mods(CONSTMODS, VARMODS)
  msfit_check_tolerance(MASS, TOL, TOLU)
  msfit_check_matches(MINMATCH, MINPMATCH, MAXMODS, MAXHITS)
  msfit_check_display(GRAPH, MOWSE, MOWSEP)
  msfit_check_instrument(INST, FORMAT)
  
  # generate peaklist
  PEAKS <- paste(PEAKS, collapse="\n")
  
  # Create HTML inputs for contant modifications (CONSTMODS)
  const_mods_html <- if (length(CONSTMODS) > 0) {
    paste0('<input type="hidden" name="const_mod" value="', CONSTMODS, '"/>', collapse = "\n")
  } else {
    ""
  }
  
  # Create HTML inputs for variable modifications (VARMODS)
  var_mods_html <- if (length(VARMODS) > 0) {
    paste0('<input type="hidden" name="mod_AA" value="', VARMODS, '"/>', collapse = "\n")
  } else {
    ""
  }
  
  # create and submit HTML
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
    
    <input type="hidden" name="database" value="', DATABASE, '">
    <input type="hidden" name="species" value="', SPECIES, '">
    
    <input type="hidden" name="dna_frame_translation" value="', FRAME, '">
    <input type="hidden" name="n_term_aa_limit" value="', AALIMIT, '">
    
    <input type="hidden" name="output_type" value="', OUTPUT, '">
    <input type="hidden" name="sort_type" value="', SORT, '">
    <input type="hidden" name="ms_report_homologous_proteins" value="', HOMOLOGUES, '">
    
    <input type="hidden" name="enzyme" value="', CLEAVAGE, '">
    <input type="hidden" name="missed_cleavages" value="', MISSEDC, '">
    
    ', const_mods_html, '
    ', var_mods_html, '
    
    <input type="hidden" name="parent_mass_convert" value="', MASS, '">
    <input type="hidden" name="ms_parent_mass_tolerance" value="', TOL, '">
    <input type="hidden" name="ms_parent_mass_tolerance_units" value="', TOLU, '">
    <input type="hidden" name="ms_parent_mass_systematic_error" value="0">
    
    <input type="hidden" name="min_matches" value="', MINMATCH, '">
    <input type="hidden" name="min_parent_ion_matches" value="', MINPMATCH, '">
    <input type="hidden" name="ms_max_modifications" value="', MAXMODS, '">
    <input type="hidden" name="ms_max_reported_hits" value="', MAXHITS, '">
    
    <input type="hidden" name="mowse_on" value="', MOWSE, '">
    <input type="hidden" name="mowse_pfactor" value="', MOWSEP, '">
    
    <input type="hidden" name="ms_prot_low_mass" value="', LOWMASS, '">
    <input type="hidden" name="ms_prot_high_mass" value="', HIGHMASS, '">
    
    <input type="hidden" name="low_pi" value="', LOWPI, '">
    <input type="hidden" name="high_pi" value="', HIGHPI, '">
    <input type="hidden" name="full_pi_range" value="', RANGEPI, '">
    
    <input type="hidden" name="instrument_name" value="', INST, '">
    <input type="hidden" name="data_format" value="', FORMAT, '">
    
    <input type="hidden" name="display_graph" value="', GRAPH, '">
    <input type="hidden" name="comment" value="">
    
    <input type="hidden" name="output_filename" value="results">
    
    <textarea name="user_protein_sequence" rows="6" cols="70">', SEQUENCE, '</textarea>
    <textarea name="data">', PEAKS, '</textarea>
    
    </form>
    </body>
    </html>'
  )
  
  temp <- tempfile(fileext = ".html")
  writeLines(text = html, con = temp)
  
  if (TEST == FALSE){
    browseURL(url = temp)
  }
}


# =================================================== #
# ===== ASSISTING FUNCTIONS ========================= #
# =================================================== #

msfit_check_peaks <- function(PEAKS){
  if (length(PEAKS) < 1){
    stop("Argument `PEAKS` cannot be empty.")
  } else if (length(PEAKS) < 2 || !is.numeric(PEAKS)){
    stop("Argument `PEAKS` needs to be a numeric vector.")
  }
}


msfit_check_sequence <- function(SEQUENCE){
  aminos <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
              "Q", "R", "S", "T", "V", "W", "Y")

  if(!is.null(SEQUENCE)){
    if (length(SEQUENCE) != 1){
      warning("Argument `SEQUENCE` should be empty or one sequence.")
    } else if (is.na(SEQUENCE)){
      warning("Argument `SEQUENCE` should be empty or one sequence.")
    } else if (!is.character(SEQUENCE)){
      warning("Argument `SEQUENCE` should be empty or one character string.")
    } else {
      chars <- unlist(strsplit(SEQUENCE, ""))
      for (aa in chars) {
        if (!(aa %in% aminos)) {
          warning("Argument `SEQUENCE` contains unexpected characters.")
          break
        }
      }
    }
  }
}


msfit_check_database <- function(DATABASE, SPECIES){
  databases <- c("User Protein", "NCBInr.2013.6.17", "NCBInr.2013.6.17.random", 
             "NCBInr.2013.6.17.random.concat", "NCBInr.refseq.2009.09.14", 
             "NCBInr.refseq.2009.09.14.random", 
             "NCBInr.refseq.2009.09.14.random.concat", 
             "NCBInr.refseq.2009.09.14.random.sub", 
             "NCBInr.refseq.2009.09.14.sub", "nextprot_all", 
             "nextprot_all.random", "nextprot_all.random.concat", 
             "Pdefault.TAIR10_pep_20101214_updated", 
             "Pdefault.TAIR10_pep_20101214_updated.random", 
             "Pdefault.TAIR10_pep_20101214_updated.random.concat", 
             "Refseq.2024.11.08.vertebrate_mammalian", 
             "Refseq.2024.11.08.vertebrate_mammalian.random", 
             "Refseq.2024.11.08.vertebrate_mammalian.random.concat", 
             "SwissProt.2016.9.6", "SwissProt.2016.9.6.random", 
             "SwissProt.2016.9.6.random.concat", "SwissProt.2017.11.01", 
             "SwissProt.2017.11.01.random", 
             "SwissProt.2017.11.01.random.concat", "SwissProt.2020.09.02", 
             "SwissProt.2020.09.02.random", 
             "SwissProt.2020.09.02.random.concat", "SwissProt.2021.06.18", 
             "SwissProt.2021.06.18.random", 
             "SwissProt.2021.06.18.random.concat", "SwissProt.2025.02.26", 
             "SwissProt.2025.02.26.random", 
             "SwissProt.2025.02.26.random.concat", 
             "SwissProt.Human.2024.02.02", "SwissProt.Human.2024.02.02.random", 
             "SwissProt.Human.2024.02.02.random.concat", 
             "UniProtKB.2017.11.01", "UniProtKB.2017.11.01.random", 
             "UniProtKB.2017.11.01.random.concat", "UniProtKB.2020.09.02", 
             "UniProtKB.2020.09.02.random", "UniProtKB.2020.09.02.random.concat")
  
  speciess <- c("All", "HUMAN MOUSE", "HUMAN RODENT", "MODEL PLANTS", 
                "RODENT", "ROACH LOCUST BEETLE", "MICROORGANISMS", 
                "ARABIDOPSIS", "ARABIDOPSIS THALIANA", 
                "ARCHAEOGLOBUS FULGIDUS", "BACILLUS SUBTILIS", "BACULOVIRIDAE", 
                "BORRELIA BURGDORFERI", "BOS TAURUS", "BOTHROPS INSULARIS", 
                "BRASSICA", "BRUGIA MALAYI", "CAENORHABDITIS ELEGANS", 
                "CANDIDA", "CANIS FAMILIARIS", "CAPRA HIRCUS", 
                "CERCOPITHECUS AETHIOPS", "CHLAMYDIA TRACHOMATIS", 
                "CHLAMYDOPHILA PNEUMONIAE", "DANIO RERIO", 
                "DICTYOSTELIUM DISCOIDEUM", "DROSOPHILA MELANOGASTER", 
                "EQUUS CABALLUS", "ESCHERICHIA COLI", "FELIS CATUS", 
                "FRANCISELLA TULARENSIS", "GALLUS GALLUS", "GLYCINE MAX", 
                "GREEN PLANTS", "GORILLA GORILLA", "HAEMOPHILUS", 
                "HAEMOPHILUS INFLUENZAE", "HELICOBACTER PYLORI", 
                "HOMO SAPIENS", "HORDEUM", "HUMAN ADENOVIRUS 5", "MACACA", 
                "MAMMALS", "METHANOCOCCUS JANNASCHII", "MUS MUSCULUS", 
                "MYCOBACTERIUM", "MYCOPLASMA", "NATRIALBA MAGADII", 
                "NEISSERIA", "ORYCTOLAGUS CUNICULUS", "ORYZA", "ORYZA SATIVA", 
                "OVIS ARIES", "PAN TROGLODYTES", "PSEUDOMONAS AERUGINOSA", 
                "RATTUS NORVEGICUS", "SACCHAROMYCES", 
                "SACCHAROMYCES CEREVISIAE", "SALMONELLA", 
                "SCHIZOSACCHAROMYCES POMBE", "SECALE", 
                "SHEWANELLA ONEIDENSIS", "SUS SCROFA", "SYNECHOCYSTIS", 
                "TOXOPLASMA GONDII", "TREPONEMA PALLIDUM", "TRITICUM", 
                "TRYPANOSOMA", "UNKNOWN", "UNREADABLE", "XENOPUS LAEVIS", "ZEA")
  
  if (length(DATABASE) != 1) {
    stop("Argument `DATABASE` should be one database as a character string.")
  } else if (!(DATABASE %in% databases)) {
    stop("Argument `DATABASE` contains unexpected database.")
  }
  
  if (length(SPECIES) != 1) {
    stop("Argument `SPECIES` should be one species as a character string.")
  } else if (!(SPECIES %in% speciess)) {
    stop("Argument `SPECIES` contains unexpected species")
    
  }
}


msfit_check_frame <- function(FRAME, AALIMIT){
  # TODO
}


msfit_check_output <- function(OUTPUT, SORT, HOMOLOGUES){
  # TODO
}


msfit_check_prerequisites <- function(LOWMASS, HIGHMASS, LOWPI, HIGHPI, RANGEPI){
  # TODO
}


msfit_check_digest <- function(CLEAVAGE, MISSEDC){
  enzymes <- c("Trypsin", "TrypsinPro", "SlymotrypsinFYWKR", "Chymotrypsin", 
               "Chymotrypsin FYW", "ChymotrypsinFWYMEDLN", "V8 DE", "V8 E", 
               "Lys-C", "Lys-C-Pro", "Lys-N", "Arg-C", "Asp-N", "Asp-C", 
               "DE-N", "CNBr", "Pepsin(porcine gastric)", "Glu-C", "Tyr-C", 
               "Thermolysin", "Elastase", "a-Lytic Protease (aLP)", 
               "Full Protein", "Pro-C", "Pro-N", "Proteinase K", 
               "Rhizopuspepsin", "Cys-N", "KR-N", "Hydroxylamine", 
               "CNBr/Trypsin", "CNBr/V8 DE", "CNBr/V8 E", "CNBr/Lys-C", 
               "CNBr/Asp-N", "CNBr/Asp-C", "CNBr/Arg-C", "CNBr/Trypsin/V8 DE", 
               "CNBr/Trypsin/V8 E", "CNBr/Arg-C/V8 E", "CNBr/Lys-C/V8 E", 
               "Trypsin/Asp-N", "Trypsin/DE-N", "Trypsin/V8 DE", 
               "Trypsin/V8 E", "Trypsin/Glu-C", "Trypsin/Chymotrypsin", 
               "Trypsin/Pepsin(porcine gastric)", "Chymotrypsin/Arg-C", 
               "Lys-C/Trypsin", "Lys-C/V8 DE", "Lys-C/V8 E", "Lys-C/Asp-N", 
               "Lys-C/DE-N", "Lys-C/Glu-C", "V8 DE/Chymotrypsin", "Arg-C/V8 E", 
               "Glu-C/Asp-N", "Glu-C/Chymotrypsin", "DE-N/Cys-N", 
               "Asp-N/Asp-C", "Asp-C/Cys-N", "Lys-C/Lys-N", "Trypsin/KR-N", 
               "Trypsin/Hydroxylamine")

  if (length(CLEAVAGE) != 1){
    stop("Argument `CLEAVAGE` should contains one enzyme.")
  }else if (!(CLEAVAGE %in% enzymes)) {
    stop("Argument `CLEAVAGE` contains unexpectet enzyme.")
  } else if (length(MISSEDC) != 1){
    warning("Argument `MISSEDC` should contains one integer.")
  } else if (!is.numeric(MISSEDC)) {
    warning("Argument `MISSEDC` should be an integer.")
  }else if (!(MISSEDC %in% c(0:9))) {
    warning("Argument `MISSEDC` should be an integer between '0' and '9'.")
  }
}


msfit_check_mods <- function(CONSTMODS, VARMODS){
  # TODO
}


msfit_check_tolerance <- function(MASS, TOL, TOLU){
  # TODO
}


msfit_check_matches <- function(MINMATCH, MINPMATCH, MAXMODS, MAXHITS){
  # TODO
}


msfit_check_display <- function(GRAPH, MOWSE, MOWSEP){
  # TODO
}


msfit_check_instrument <- function(INST, FORMAT){
  # TODO
}
