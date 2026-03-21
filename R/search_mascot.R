# TODO: go through website HTML and orter arguments and local html

# =================================================== #
# ===== MAIN FUNCTION =============================== #
# =================================================== #

search_mascot <- function(
    PEAKS,
  
    USERNAME,
    USEREMAIL,
    
    SEARCH = "PMF",
    FORMVER = 1.01,
    
    DATABASE = "SwissProt",
    
    CLEAVAGE = "Trypsin", 
    MISSEDC = 1,
    
    CONSTMODS = NULL,
    VARMODS = "Oxidation (M)",
    
    MASS = "Monoisotopic",
    CHARGE = "1+",
    TOL = 0.3,
    TOLU = "Da",
    
    REPORT = "AUTO",
    
    TEST = FALSE
  ){
  
  #check inputs
  mascot_check_peaks(PEAKS)
  mascot_check_auth(USERNAME, USEREMAIL)
  mascot_check_database(DATABASE)
  mascot_check_digest(CLEAVAGE, MISSEDC)
  mascot_check_mods(CONSTMODS, VARMODS)
  mascot_check_tolerance(MASS, CHARGE, TOL, TOLU)
  mascot_check_display(REPORT)

  # generate peaklist
  PEAKS <- paste(PEAKS, collapse="\n")
  
  # Create HTML inputs for contant modifications (CONSTMODS)
  const_mods_html <- if (length(CONSTMODS) > 0) {
    paste0('<input type="hidden" name="MODS" value="', CONSTMODS, '"/>', collapse = "\n")
  } else {
    ""
  }
  
  # Create HTML inputs for variable modifications (VARMODS)
  var_mods_html <- if (length(VARMODS) > 0) {
    paste0('<input type="hidden" name="IT_MODS" value="', VARMODS, '"/>', collapse = "\n")
  } else {
    ""
  }
  
  # create and submit HTML
  html <- paste0(
    
    '<!DOCTYPE html>
    <html>
    <body onload="document.forms[0].submit()">
    
    <form method="POST" action="https://www.matrixscience.com/cgi/nph-mascot.exe?1+-batch+-format+msr" enctype="multipart/form-data">
    
    <input type="hidden" name="USERNAME" value="', USERNAME, '">
    <input type="hidden" name="USEREMAIL" value="', USEREMAIL, '">
    
    <input type="hidden" name="SEARCH" value="', SEARCH, '">
    <input type="hidden" name="FORMVER" value="', FORMVER, '">
    
    <input type="hidden" name="DB" value="', DATABASE, '">
    
    <input type="hidden" name="CLE" value="', CLEAVAGE, '">
    <input type="hidden" name="PFA" value="', MISSEDC, '">
    
    ', const_mods_html, '
    ', var_mods_html, '
    
    <input type="hidden" name="TOL" value="', TOL, '">
    <input type="hidden" name="TOLU" value="', TOLU, '">
    
    <input type="hidden" name="MASS" value="', MASS, '">
    <input type="hidden" name="CHARGE" value="', CHARGE, '">
    
    <input type="hidden" name="REPORT" value="', REPORT, '">
    
    <textarea name="QUE">', PEAKS, '</textarea>
    
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

mascot_check_peaks <- function(PEAKS){
  if (length(PEAKS) < 1){
    stop("Argument `PEAKS` cannot be empty.")
  } else if (length(PEAKS) < 2 || !is.numeric(PEAKS)){
    stop("Argument `PEAKS` needs to be a numeric vector.")
  }
}


mascot_check_auth <- function(USERNAME, USEREMAIL){
  if (length(USERNAME) < 1){
    stop("Argument `USERNAME` cannot be empty.")
  } else if (!is.character(USERNAME)){
    stop("Argument `USERNAME` can only contain characters.")
  } else if (!grepl("^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.\\-]+\\.[a-zA-Z]{2,}$", USEREMAIL)){
    stop("Argument `USEREMAIL` needs to be a valid mail address.")
  }
}


mascot_check_database <- function(DATABASE){
  databases <- c("Fungi_EST", "Human_EST", "Invertebrates_EST", "Mammals_EST", 
                 "Mus_EST", "Plants_EST", "Prokaryotes_EST", "Rodents_EST", 
                 "Vertebrates_EST", "contaminants", "cRAP", "SwissProt", 
                 "UP186698_X_laevis", "UP19116_T_aestivum", "UP1940_C_elegans", 
                 "UP2195_D_discoideum", "UP219602_F_oxysporum", 
                 "UP2311_S_cerevisiae", "UP241690_T_harzianum", 
                 "UP2485_S_pombe", "UP2494_R_norvegicus", "UP254291_M_gilvum", 
                 "UP279841_T_thermophilus", "UP317484_G_aquaeductus", 
                 "UP437_D_rerio", "UP5226_T_rubripes", "UP5640_H_sapiens", 
                 "UP589_M_musculus", "UP59680_O_sativa", "UP625_E_coli_K12", 
                 "UP6548_A_thaliana", "UP6906_C_reinhardtii", "UP7305_Z_mays", 
                 "UP803_D_melanogaster", "UP808_M_pneumoniae", 
                 "UP8227_S_scrofa", "UP8311_R_communis", "UP9136_B_taurus")

  if (length(DATABASE) != 1){
    stop("Argument `DATABASE` should contains one database.")
  }else if (!(DATABASE %in% databases)) {
    stop("Argument `DATABASE` contains unexpectet database.")
  }
}


mascot_check_digest <- function(CLEAVAGE, MISSEDC){
  enzymes <- c("Trypsin", "Trypsin/P", "Arg-C", "Asp-N", "Asp-N_ambic", 
               "Chymotrypsin", "CNBr", "CNBr+Trypsin", "Formic_acid", "Lys-C", 
               "Lys-C/P", "LysC+AspN", "Lys-N", "PepsinA", "semiTrypsin", 
               "TrypChymo", "TrypsinMSIPI", "TrypsinMSIPI/P", "V8-DE", "V8-E")

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


mascot_check_mods <- function(CONSTMODS, VARMODS){
  # TODO
}


mascot_check_tolerance <- function(MASS, CHARGE, TOL, TOLU){
  if (MASS != "Monoisotopic"){
    stop("Argument `MASS` should be 'Monoisotopic'.")
  } else if (CHARGE != "1+") {
    stop("Argument `CHARGE` should be '1+'.")
  } else if (length(TOL) != 1){
    stop("Argument `TOL` should contains one floating point value.")
  } else if (!is.numeric(TOL)) {
    stop("Argument `TOL` should be a floating point value.")
  } else if (TOL < 0) {
    stop("Argument `TOL` should be a positive floating point value.")
  } else if (length(TOLU) != 1){
    stop("Argument `TOLU` should contains one unit.")
  } else if (!(TOLU %in% c("Da", "ppm", "mmu", "%"))) {
    stop("Argument `TOLU` should be an integer.")
  }
}


mascot_check_display <- function(REPORT){
  # TODO
}
