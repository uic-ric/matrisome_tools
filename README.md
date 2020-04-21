# Matrisome Project tools

The following are tools used by the Matrisome project (http://matrisome.org) to process proteomic data to add Matrisome annotations to proteins and generate a summary report.

## process_ecm.sh

Tool will process mzML and idXML/pepXML file to obtain PSM counts and EIC intensities for each peptide then annotate proteins using ECM database

```console
usage: process_ecm.sh [mzML file] [idXML/pepXML file] [ECM DB file] [output peptide table] [output protein table]
```

## ecm_report.R

Script will generate a summary report of output protein tables from the process_ecm.sh script.  The multiple protein tables from different samples can be provided to create an overall report of a given sample set.  The script will generate the following...

- **Summary report** - Summary of number of proteins, number of peptides, number of PSMs/spectra, and EIC sum for each sample broken down for each Matrisome division as well as collagens.  This output file is specified with the `-o/--output` option.
- **Combined sample report** - Report of each protein with Matrisome annotations and summed counts of peptides, PSMs/spectra and EIC area for each sample.
- **Plots** - Pie or barcharts of summariezed protein, peptide, PSM/spectra counts and EIC area colored by Matrisome division.

```console


Usage: /export/home/clustcrilab/internal_tools/george_tools/matrisome_qc/ecm_report.R [options]
Generate report for ECM data.

Options:
	-o OUTPUT, --output=OUTPUT
		Output summary report (text). Default is STDOUT.

	-d DETAILS, --details=DETAILS
		Output combined protein report (text).

	--gene_details=GENE_DETAILS
		Output combined protein report summarized by gene (text).

	-i MANIFEST, --manifest=MANIFEST
		Manifest of input files. Assumes a header row. If just a list of files, will assume there is NOT a header. If this is not specified then provide files as arguments.

	-m MIN_PEPTIDES, --min_peptides=MIN_PEPTIDES
		Minimum number of peptides required to include a protein. Default is 2

	-u, --unique_peptides
		Minimum peptides must be unique.

	-p PLOT, --plot=PLOT
		Output plot file(s).

	-t TYPE, --type=TYPE
		Plot format, either png or pdf. Default is pdf.

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

## Galaxy tools

### Building Docker images

The current Galaxy tool XML files are written to use Docker containers.  The following Dockerfiles are provided for each of the Galaxy tools

* Dockerfile.process_ecm
* Dockerfile.ecm_report

To build the Docker images execute the following, where `tool` is process_ecm or ecm_report.

```console
docker build -f Dockerfile.<tool> -t <tool> ./
```

*Note*: The process_ecm.py script fetches protein information from Uniprot (https://uniprot.org).  When running a `process_ecm` Docker container, the container needs to have Internet access.  If you are using the `local_docker` job runner be sure to set the `docker_net` parameter for the `local_docker` destination.  The following is an example enabling network access via a bridge interface.

```xml
<destination id="local_docker" runner="local">
   ...
   <param id="docker_net">bridge</param>
   ...
</destination>
```



