# MassSpectFPAutoTool

Automated preprocessing and peptide mass fingerprint (PMF) analysis for MALDI-TOF mass spectrometry data in **R**.

The project contains two scripts that process MALDI spectra, detect monoisotopic peaks, remove common contaminants, and optionally submit peak lists to protein identification services.

---

## Scripts

* **MassSpectFPAutoTool for Mascot.R**
  Processes spectra and submits peak lists to **Mascot**.

* **MassSpectFPAutoTool for MS-Fit.R**
  Processes spectra and submits peak lists to **MS-Fit (Protein Prospector)**.
  Optionally performs an additional `.random.concat` search for FDR estimation.

---

## Requirements

* R (recommended ≥ 4.0)
* Internet connection for database submission

Required R packages (installed automatically if missing):

* `BiocManager`
* `MALDIquant`
* `MALDIquantForeign`

---

## Workflow

Both scripts perform the following steps:

1. Import MALDI spectrum
2. Preprocess spectrum
3. Detect peaks
4. Remove trypsin and keratin contaminant peaks
5. Extract monoisotopic peaks
6. Automatically adjust SNR to obtain ~20–40 peaks
7. Export results
8. Optionally submit search to Mascot or MS-Fit

---

## How to Run

1. Open the script in **R** or **RStudio**.

    - for easy file selection, the file should be in the same directory as the spectrum files.

2. Run the script.

    - In **RStudio**:
        - Select all code (`Cmd/Ctrl + a`)
        - Run the script (`Cmd/Ctrl + Shift + Enter`)

3. Select a spectrum file when prompted.

    - for **Mascot** search, only one file can be selected as the free trial server does not allow batch processing of files. **Please adhere to this restriction even when editing the code.**

4. Choose whether to run the database search.

    - for **Mascot** you need to specify your USERNAME = Name and USEREMAIL = email. The first time you will be sent a mail to verify your mail address.

The script will then process the spectrum(s) automatically.

---

## Output

A folder named **`results`** is created in the input file directory.

Generated files include:

* Original spectrum plot
* Preprocessed spectrum plot
* Detected peak plot
* Monoisotopic peak plot
* Monoisotopic peak list (`.txt`)
* Processing log (`processing_log.txt`)

---

## Parameters

Search parameters can be modified directly inside the scripts:

* `mascot_parameters` for Mascot
* `msfit_parameters` for MS-Fit

Adjust these values to match experimental conditions.

---

## Notes

* Designed for MALDI-TOF peptide mass fingerprinting.
* Optimal results are obtained with spectra producing **20–40 monoisotopic peaks**.
* Contaminant peaks from trypsin autolysis and keratin are automatically removed.

---

## Disclaimer

This software is provided **for educational and research purposes only**.

* No guarantee of correctness or suitability
* Not intended for clinical, diagnostic, or commercial use
* Predictions should always be validated independently

Use at your own risk.

---

## Author

Spike Murphy Müller

---

## Copyright

© 2026 Spike Murphy Müller · Licensed under the MIT License
