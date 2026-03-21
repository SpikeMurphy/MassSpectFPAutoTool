\name{search_msfit}
\alias{search_msfit}
\title{Perform a Protein Prospector MS-Fit Peptide Mass Fingerprint Search from a Peaklist}
\description{This function takes a peak list as input and submits a request 
to a Protein Prospector MS-Fit server to perform a peptide mass fingerprint 
(PMF) database search using specified parameters. Results are opened in the 
user's default browser.

\strong{Important:} MS-Fit servers do not allow batch submissions. This function 
should only be used for individual peak lists and must not be used to automate 
the analysis of multiple peak lists.
}
\usage{
search_msfit(
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
    
    CLEAVAGE = "Trypsin",
    MISSEDC = 1,
    
    CONSTMODS = NULL,
    VARMODS = "Oxidation of M",
    
    MASS = "monoisotopic",
    TOL = 0.3,
    TOLU = "Da",
    
    MINMATCH = 4,
    MINPMATCH = 1,
    MAXMODS = 1,
    MAXHITS = 5,
    
    GRAPH = 1,
    MOWSE = 1,
    MOWSEP = 0.4,
    
    INST = "MALDI-TOFTOF",
    FORMAT = "PP M/Z Charge"
)
}
\arguments{
  \item{PEAKS}{Numeric vector. Mass-to-charge (m/z) values representing the peak list used for the peptide mass fingerprint search.}
  
  \item{SEQUENCE}{Character string. Optional protein or nucleotide sequence to restrict the search. The default argument is \code{"NULL"}.}
  
  \item{DATABASE}{Character string. Character string. Protein sequence database to search against. The default argument is \code{"SwissProt.2025.02.26"}.}
  
  \item{SPECIES}{Character string. Taxonomic restriction for the search. The default argument is \code{"All"}.}
  
  \item{FRAME}{Integer. Reading frame for DNA translation. The default argument is \code{3}.}
  
  \item{AALIMIT}{Character string. Amino terminal amino acid limit. The default argument is \code{"NULL"}.}

  \item{OUTPUT}{Character string. Output format. The default argument is \code{"HTML"} with an alternative being \code{"XML"}.}
  
  \item{SORT}{Character string. Sorting method for results. The default argument is \code{"Score Sort"}.}
  
  \item{HOMOLOGUES}{Character string. Homologues to report.The default argument is \code{"Interesting"} with alternatives being \code{"All"} or \code{"None"}.}
  
  \item{LOWMASS}{Numeric. Lower molecular weight limit for proteins.The default argument is \code{1000}.}
  
  \item{HIGHMASS}{Numeric. Upper molecular weight limit for proteins. The default argument is \code{125000}.}
  
  \item{LOWPI}{Numeric. Lower isoelectric point (pI) limit. The default argument is \code{3}.}
  
  \item{HIGHPI}{Numeric. Upper isoelectric point (pI) limit. The default argument is \code{10}.}
  
  \item{RANGEPI}{Character string. Flag to enable pI range filtering. The default argument is \code{"1"}.}
  
  \item{CLEAVAGE}{Character string. Enzyme used for protein digestion. The default argument is \code{"Trypsin"}.}
  
  \item{MISSEDC}{Integer. Number of allowed missed cleavages. The default argument is \code{1}.}
  
  \item{CONSTMODS}{Character string. Fixed modifications. The default argument is \code{"NULL"}.}
  
  \item{VARMODS}{Character string. Variable modifications. The default argument is \code{"Oxidation of M"}.}

  \item{MASS}{Character string. Mass type used. The default argument is \code{"monoisotopic"} referring to the first peak of an isotopic distribution.}
  
  \item{TOL}{Numeric. Peptide mass tolerance value. The default argument is \code{0.3}.}
  
  \item{TOLU}{Character string. Units for mass tolerance. The default argument is \code{"Da"} for Dalton with an alternative being \code{"ppm"} for parts per million.}

  \item{MINMATCH}{Integer. Minimum number of peptide matches required. The default argument is \code{4}.}
  
  \item{MINPMATCH}{Integer. Minimum number of parent mass matches required. The default argument is \code{1}.}
  
  \item{MAXMODS}{Integer. Maximum number of modifications per peptide. The default argument is \code{1}.}
  
  \item{MAXHITS}{Integer. Maximum number of hits returned. The default argument is \code{5}.}
  
  \item{GRAPH}{Integer (0/1). Whether to display graphical output. Chose between \code{0} for FALSE or \code{1} for TRUE. The default argument is \code{1}.}
  
  \item{MOWSE}{Integer (0/1). Whether to calculate MOWSE scores. Chose between \code{0} for FALSE or \code{1} for TRUE. The default argument is \code{1}.}
  
  \item{MOWSEP}{Numeric. MOWSE probability threshold. The default argument is \code{0.4}.}
  
  \item{INST}{Character string. Instrument type. The default argument is \code{"MALDI-TOFTOF"}.}
  
  \item{FORMAT}{Character string. Input peak list format. The default argument is \code{"PP M/Z Charge"}.}
}
\value{
No return value. As a side effect, the function opens the user's default browser 
and submits an HTML-form to the MS-Fit server to initiate the search.
}
\examples{
\dontrun{
search_msfit(
  PEAKS = c(
      # human serum albumin
      802.4781, 809.4992, 817.4778, 854.4519, 875.5098, 876.4971, 880.4411, 890.5128,
      892.4371, 915.4175, 927.4934, 940.4483, 951.4418, 960.5625, 973.5214, 984.4884,
      1000.6037, 1002.5578, 1013.4244, 1013.599, 1017.5364, 1018.4775, 1019.5785,
      1024.455, 1055.5884, 1074.5426, 1083.5946, 1128.6987, 1138.5674, 1140.5354,
      1141.5194, 1141.6939, 1149.5759, 1149.615, 1157.6538, 1191.5748, 1226.6052,
      1252.6572, 1257.5238, 1296.7046, 1311.7419, 1320.4905, 1336.4854, 1341.7559,
      1342.6348, 1352.6661, 1352.7685, 1358.6297, 1381.5334, 1384.5355, 1386.6206,
      1432.7538, 1442.6104, 1467.8431, 1490.9053, 1511.8428, 1548.6821, 1600.7312,
      1623.7876, 1623.9581, 1627.7996, 1639.7825, 1639.9377, 1650.8949, 1657.7751,
      1676.7481, 1692.743, 1742.8941, 1779.8887, 1795.8836, 1840.915, 1853.9102,
      1856.9099, 1875.0157, 1884.9008, 1898.9952, 1898.9952, 1915.7731, 1924.0862,
      1939.9079, 1953.9231, 1959.7881, 1975.783, 2037.0844, 2045.0954, 2177.9698,
      2193.9647, 2203.0012, 2296.2199, 2300.1056, 2312.2148, 2316.1005, 2347.0039,
      2383.0541, 2404.1709, 2414.0533, 2420.1659, 2433.2635, 2488.1449, 2496.3438,
      2515.1327, 2542.2759, 2560.272, 2576.267, 2593.2425, 2609.2374, 2686.2269,
      2721.3375, 2737.3324, 2859.3474, 2917.323, 2933.3179, 3059.4999, 3075.4948,
      3362.5225, 3365.6875, 3378.5174, 3407.6109, 3514.673, 3530.6679, 3563.8607,
      3626.8087
    )
)

search_msfit(
  PEAKS = c(
      # human serum albumin
      802.4781, 809.4992, 817.4778, 854.4519, 875.5098, 876.4971, 880.4411, 890.5128,
      892.4371, 915.4175, 927.4934, 940.4483, 951.4418, 960.5625, 973.5214, 984.4884,
      1000.6037, 1002.5578, 1013.4244, 1013.599, 1017.5364, 1018.4775, 1019.5785,
      1024.455, 1055.5884, 1074.5426, 1083.5946, 1128.6987, 1138.5674, 1140.5354,
      1141.5194, 1141.6939, 1149.5759, 1149.615, 1157.6538, 1191.5748, 1226.6052,
      1252.6572, 1257.5238, 1296.7046, 1311.7419, 1320.4905, 1336.4854, 1341.7559,
      1342.6348, 1352.6661, 1352.7685, 1358.6297, 1381.5334, 1384.5355, 1386.6206,
      1432.7538, 1442.6104, 1467.8431, 1490.9053, 1511.8428, 1548.6821, 1600.7312,
      1623.7876, 1623.9581, 1627.7996, 1639.7825, 1639.9377, 1650.8949, 1657.7751,
      1676.7481, 1692.743, 1742.8941, 1779.8887, 1795.8836, 1840.915, 1853.9102,
      1856.9099, 1875.0157, 1884.9008, 1898.9952, 1898.9952, 1915.7731, 1924.0862,
      1939.9079, 1953.9231, 1959.7881, 1975.783, 2037.0844, 2045.0954, 2177.9698,
      2193.9647, 2203.0012, 2296.2199, 2300.1056, 2312.2148, 2316.1005, 2347.0039,
      2383.0541, 2404.1709, 2414.0533, 2420.1659, 2433.2635, 2488.1449, 2496.3438,
      2515.1327, 2542.2759, 2560.272, 2576.267, 2593.2425, 2609.2374, 2686.2269,
      2721.3375, 2737.3324, 2859.3474, 2917.323, 2933.3179, 3059.4999, 3075.4948,
      3362.5225, 3365.6875, 3378.5174, 3407.6109, 3514.673, 3530.6679, 3563.8607,
      3626.8087
    ),
    SEQUENCE = "",
    
    DATABASE = "SwissProt.2025.02.26",
    SPECIES = "All",
    
    FRAME = 3,
    AALIMIT = "",

    OUTPUT = "HTML",
    SORT = "Score Sort",
    HOMOLOGUES = "Interesting",
    
    LOWMASS = 1000,
    HIGHMASS = 125000,
    LOWPI = 3.0,
    HIGHPI = 10.0,
    RANGEPI = "1",
    
    CLEAVAGE = "Trypsin",
    MISSEDC = 1,
    
    CONSTMODS = "",
    VARMODS = "Oxidation of M",
    
    MASS = "monoisotopic",
    TOL = 0.3,
    TOLU = "Da",
    
    MINMATCH = 4,
    MINPMATCH = 1,
    MAXMODS = 1,
    MAXHITS = 5,
    
    GRAPH = 1,
    MOWSE = 1,
    MOWSEP = 0.4,
    
    INST = "MALDI-TOFTOF",
    FORMAT = "PP M/Z Charge"
)
}
}