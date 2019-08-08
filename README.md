# Matrisome Project tools

The following are tools used by the Matrisome project (http://matrisome.org) to process proteomic data to add Matrisome annotations to proteins and generate a summary report.

## process_ecm.sh

Tool will process mzML and idXML/pepXML file to obtain PSM counts and EIC intensities for each peptide then annotate proteins using ECM database

```bash
usage: process_ecm.sh [mzML file] [idXML/pepXML file] [ECM DB file] [output peptide table] [output protein table]
```

## ecm_report.R

Script will generate a summary report of output protein tables from the process_ecm.sh script.  The multiple protein tables from different samples can be provided to create an overall report of a given sample set.  The script will generate the following...

- **Summary report** - Summary of number of proteins, number of peptides, number of PSMs/spectra, and EIC sum for each sample broken down for each Matrisome division as well as collagens.  This output file is specified with the `-o/--output` option.
- **Combined sample report** - Report of each protein with Matrisome annotations and summed counts of peptides, PSMs/spectra and EIC area for each sample.
- **Plots** - Pie or barcharts of summariezed protein, peptide, PSM/spectra counts and EIC area colored by Matrisome division.

```bash
Usage: ecm_report.R [options] [protein tables...]
Generate report for ECM data.

Options:
        -o OUTPUT, --output=OUTPUT
                Output summary report (text). Default is STDOUT.

        -d DETAILS, --details=DETAILS
                Output combined sample report (text).

        -m MIN_PEPTIDES, --min_peptides=MIN_PEPTIDES
                Minimum number of peptides. Default is 2

        -p PLOT, --plot=PLOT
                Output plot file.

        -t TYPE, --type=TYPE
                Output type, either png or pdf. Default is pdf.

        -s STYLE, --style=STYLE
                Plot style, either bar or pie. Default is bar.

        -h NUMBER, --height=NUMBER
                Height of plots (600px if PNG, 7in if PDF).

        -w NUMBER, --width=NUMBER
                Width of plots (600px if PNG, 7in if PDF).

        -r, --rotate
                Rotate plots (x on the vertical).

        --help
                Display this help message and exit.
```
