#!/usr/bin/env Rscript
################################################################################
# Script : ecm_report.R
# Author : George Chlipala
# Created: Feb 13, 2019
# -- Description ----------------------------------------
# Script to process output from ecm_process.py to generate a report
# -- Requirements ---------------------------------------
# R
################################################################################

library(optparse)
library(reshape2)

##
# Setup script arguments
option_list <- list( 
    make_option(c("-o", "--output"), help="Output summary report (text). Default is STDOUT.", default="-"),
    make_option(c("-d", "--details"), help="Output combined protein report (text)."),
    make_option(c("-i", "--manifest"), help=paste("Manifest of input files. Assumes a header row. If just a list of files, will assume there is NOT a header.", 
	"If this is not specified then provide files as arguments.")),
    make_option(c("-m", '--min_peptides'), help="Minimum number of peptides required to include a protein. Default is %default", type="integer", default=2),
    make_option(c("-u", "--unique_peptides"), help="Minimum peptides must be unique.", default=FALSE, action="store_true"),
    make_option(c("-p", "--plot"), help="Output plot file(s)."),
    make_option(c("-t", "--type"), help="Plot format, either png or pdf. Default is %default.", default="pdf"),
    make_option(c("-s", "--style"), help="Plot style, either bar or pie. Default is %default.", default="bar"),
    make_option(c("-h", "--height"), type="double", 
	help="Height of plots (600px if PNG, 7in if PDF).", metavar="number"),
    make_option(c("-w", "--width"), type="double", 
	help="Width of plots (600px if PNG, 7in if PDF).", metavar="number"),
    make_option(c("-r", "--rotate"), help="Rotate plots (x on the vertical).", 
	default=FALSE, action="store_true"),
    make_option(c("-?", "--help"), help="Display this help message and exit.", 
	default=FALSE, action="store_true")
    )

## Setup method to print usage/help statement
printHelp <- function(opts, return.code=0) {
  print_help(opts)
  q(status=return.code)
}

##
# Define default color palette
palette_text <- "name	color
Core matrisome	#4C4CFC
ECM glycoproteins	#C47CFD
Collagens	#3DA8FB
Proteoglycans	#77F9F1
Matrisome-associated	#EE7833
ECM-affiliated proteins	#F3A343
ECM regulators	#FADDBE
Secreted factors	#F3BEDD
Other	#D3D3D3
None	#D3D3D3"

# Parse the text to create a data frame
matrisome_palette_df <- read.delim(text=palette_text, as.is=T)

# convert the color column to a named vector
matrisome_palette <- matrisome_palette_df$color
names(matrisome_palette) <- matrisome_palette_df$name

# Generate argument parser
opt_parser <- OptionParser(option_list=option_list, 
	add_help_option=FALSE, description="Generate report for ECM data.")

# Parse arguments
args <- parse_args(opt_parser, print_help_and_exit = FALSE,
	positional_arguments = c(0,Inf) )

# Get options
opts <- args$options

# If help was specified, print usage and quit
if ( opts$help ) { printHelp(opt_parser) }

# Parse input files
data_files = character()
sample_ids <- NULL

if ( ! is.null(opts$manifest) ) {
	manifest <- read.delim(opts$manifest, as.is=T)
	if ( ncol(manifest) > 1 ) { 
		sample_ids <- manifest[,1]
		data_files <- manifest[,2]
	} else {
		data_files <- readLines(opts$manifest)
	}
} else {
	if ( length(args$args) == 0 ) { 
		write("ERROR: Must specify a file manifest or files", stderr())	
		printHelp(opt_parser, 1)
	}
	for ( arg in args$args ) { 
		files = Sys.glob(arg)
		if ( length(files) > 0 ) { 
			data_files = c(data_files, files)
		} else {
			data_files = c(data_files, arg)
		}
	}
}

##
# function to parse each sample file
# Arguments are..
# 1. data frame, of the input data
# 2. list of categories to include. Default is Collagens
### Columns of input data 
### ProteinID, Gene symbol, Gene name, Division, Category, 
###    Unique peptides, Non-unique peptides, Interspecies peptides, 
###    Unique spectra, Non-unique spectra, Interspecies spectra, 
###    Unique EIC, Non-unique EIC, Interspecies EIC
sampleReport <- function(x, categories=c("Collagens")) { 

	if ( "Species" %in% colnames(x) ) { 
		x$Species <- factor(x$Species)
		species_list <- levels(x$Species)
		species_data <- x[! is.na(x$Species), ]
	} else {
		species_list <- character()
	}

	process_spectra <- "Unique spectra" %in% colnames(x)
	process_eic <- "Unique EIC" %in% colnames(x)
	
	# Generate data frame to contain report
	group_names <- character()
	spec_names <- character()
	protein_count <- numeric()
	peptide_count <- numeric()
	spectra_count <- numeric()
	eic_sums <- numeric()
	
	row_index = 1

	# First, generate the report for the major divisions
	for ( division in levels(x$Division) ) { 
		for ( species in species_list ) { 
			this_species_data <- species_data[species_data$Species == species & species_data$Division == division, ]
			group_names[row_index] <- division
			spec_names[row_index] <- species
			protein_count[row_index] <- nrow(this_species_data)
			if ( protein_count[row_index] > 0 ) { 
				peptide_count[row_index] <- sum(this_species_data[, c("Unique peptides", "Non-unique peptides")])
				spectra_count[row_index] <- ifelse(process_spectra, sum(this_species_data[, c("Unique spectra", "Non-unique spectra")]), 0)
				eic_sums[row_index] <- ifelse(process_eic, sum(this_species_data[, c("Unique EIC", "Non-unique EIC")]), 0)
			} else {
				peptide_count[row_index] <- 0
				spectra_count[row_index] <- 0
				eic_sums[row_index] <- 0
			}
			row_index = row_index + 1
		}
		this_division_data <- x[x$Division == division, ]
		group_names[row_index] <- division
		spec_names[row_index] <- "ALL"
		protein_count[row_index] <- nrow(this_division_data)
		if ( protein_count[row_index] > 0 ) { 
			peptide_count[row_index] <- sum(this_division_data$Total_peptides)
			spectra_count[row_index] <- ifelse(process_spectra, sum(this_division_data$Total_spectra), 0)
			eic_sums[row_index] <- ifelse(process_eic, sum(this_division_data$Total_EIC), 0)
		} else {
			peptide_count[row_index] <- 0
			spectra_count[row_index] <- 0
			eic_sums[row_index] <- 0
		}	
		row_index = row_index + 1
	}	

	# Finally, generate the report for the selected categories
	for ( c in categories ) { 
		for ( species in species_list ) { 
			this_species_data <- species_data[species_data$Species == species & species_data$Category == c, ]
			group_names[row_index] <- c
			spec_names[row_index] <- species
			protein_count[row_index] <- nrow(this_species_data)
			if ( protein_count[row_index] > 0 ) { 
				peptide_count[row_index] <- sum(this_species_data[, c("Unique peptides", "Non-unique peptides")])
				spectra_count[row_index] <- ifelse(process_spectra, sum(this_species_data[, c("Unique spectra", "Non-unique spectra")]), 0)
				eic_sums[row_index] <- ifelse(process_eic, sum(this_species_data[, c("Unique EIC", "Non-unique EIC")]), 0)
			} else {
				peptide_count[row_index] <- 0
				spectra_count[row_index] <- 0
				eic_sums[row_index] <- 0
			}
			row_index = row_index + 1
		}

		this_category_data <- x[x$Category == c, ]
		group_names[row_index] <- c
		spec_names[row_index] <- "ALL"
		protein_count[row_index] <- nrow(this_category_data)
		if ( protein_count[row_index] > 0 ) { 
			peptide_count[row_index] <- sum(this_category_data$Total_peptides)
			spectra_count[row_index] <- ifelse(process_spectra, sum(this_category_data$Total_spectra), 0)
			eic_sums[row_index] <- ifelse(process_eic, sum(this_category_data$Total_EIC), 0)
		} else {
			peptide_count[row_index] <- 0
			spectra_count[row_index] <- 0
			eic_sums[row_index] <- 0
		}	
		row_index = row_index + 1
	}

	# Generate the totals
	for ( species in species_list ) { 
		group_names[row_index] <- "TOTAL"	
		spec_names[row_index] <- species
		this_species_data <- species_data[species_data$Species == species, ]
		protein_count[row_index] <- nrow(this_species_data)
		peptide_count[row_index] <- sum(this_species_data[, c("Unique peptides", "Non-unique peptides")])
		spectra_count[row_index] <- ifelse(process_spectra, sum(this_species_data[, c("Unique spectra", "Non-unique spectra")]), 0)
		eic_sums[row_index] <- ifelse(process_eic, sum(this_species_data[, c("Unique EIC", "Non-unique EIC")]), 0)
		row_index = row_index + 1
	}

	group_names[row_index] <- "TOTAL"	
	spec_names[row_index] <- "ALL"
	protein_count[row_index] <- nrow(x)
	peptide_count[row_index] <- sum(x$Total_peptides)
	spectra_count[row_index] <- ifelse(process_spectra, sum(x$Total_spectra), 0)
	eic_sums[row_index] <- ifelse(process_eic, sum(x$Total_EIC), 0)

	# Generate the report data frame
	report <- data.frame(Group=group_names, Species=spec_names,
		Protein_count=protein_count, Percent_protein=protein_count / protein_count[row_index],
		Peptide_count=peptide_count, Percent_peptides=peptide_count / peptide_count[row_index],
		Spectra_count=spectra_count, Percent_spectra=ifelse(process_spectra, spectra_count / spectra_count[row_index], spectra_count),
		EIC_intensity=eic_sums, Percent_EIC=ifelse(process_eic, eic_sums / eic_sums[row_index], eic_sums))

	# Return the report
	return(report)
}

# Parse over the different samples and generate the full report
full_report = NULL
sample_details = NULL
counts_data = data.frame(Sample=character(), Division=character(), ProteinID=character(), Peptides=numeric(), Spectra=numeric())

for ( s in 1:length(data_files) ) { 
	sample_file <- data_files[s]
	# Read inputfile
	ecm_data <- read.delim(sample_file, comment.char='#', check.names=F)
	ecm_data[,'Gene name'] <- as.character(ecm_data[, 'Gene name'])
	ecm_data[,'Gene symbol'] <- as.character(ecm_data[, 'Gene symbol'])
	if ( "Species" %in% colnames(ecm_data) ) { 
		ecm_data$Species <- gsub("^([A-Za-z]+ [a-z]+) .+$", "\\1", ecm_data$Species)
	}
	# If the sample_ids vector is set (a manifest was used) then use the sample name in the manifest
	if ( ! is.null(sample_ids) ) { 
		sample_name <- sample_ids[s]
	} else {
		# Check if the sample name is contained in the file header. 
		# Used in Galaxy pipeline
		# Read the comment lines, if they exist
		ecm_comments <- grep("^##", readLines(sample_file, n=10), value=T)
		if ( length(ecm_comments) > 0 ) { 
			for ( line in ecm_comments ) { 	
				parts <- unlist(strsplit(line, ':'))		
				field <- tolower(trimws(sub("^##+ *", "", parts[1])))
				value <- trimws(paste(parts[-1], collapse=':'))
				if ( field == "sample" ) { sample_name = value }
			}
		} else {
			sample_name <- sub("\\.(txt|tsv|csv|tab)$", "", basename(sample_file))
		}
	}

	ecm_data$Total_peptides <- apply(ecm_data[, c("Unique peptides", "Non-unique peptides", "Interspecies peptides")], 1, sum)

	# Remove any proteins that do not have the minimum number of peptides
	if ( opts$unique_peptides ) {  
		ecm_data <- ecm_data[ecm_data[, "Unique peptides"] >= opts$min_peptides, ]
	} else {
		ecm_data <- ecm_data[ecm_data$Total_peptides >= opts$min_peptides, ]
	}

	# Get counts for each protein.  For histograms in report
	sample_counts = data.frame(Sample=rep(sample_name, nrow(ecm_data)), 
		ecm_data[,c("ProteinID", "Division")], Peptides=ecm_data$Total_peptides)

	# Compute total spectra counts (if included in the data)
	if ( "Unique spectra" %in% colnames(ecm_data) ) { 
		ecm_data$Total_spectra <- apply(ecm_data[, c("Unique spectra", "Non-unique spectra", "Interspecies spectra")], 1, sum)		
		sample_counts$Spectra <- ecm_data$Total_spectra
	} else {
		sample_counts$Spectra <- NA
	}

	# Compute total EIC counts (if included in the data)
	if ( "Unique EIC" %in% colnames(ecm_data) ) { 
		ecm_data$Total_EIC <- apply(ecm_data[, c("Unique EIC", "Non-unique EIC", "Interspecies EIC")], 1, sum)		
	}

	# Add basic counts information for this sample to the counts_data data frame
	counts_data = rbind(counts_data, sample_counts)

	# Generate the report for this sample
	a_report <- sampleReport(ecm_data)
	# Melt the report table
	a_report <- melt(a_report, id.vars=c("Group", "Species"))
	a_report$column <- paste(a_report$variable, a_report$Species, a_report$Group, sep=":")		
	# Set the value column to the sample name
	colnames(a_report)[4] <- sample_name
	# Add this sample to the full report
	print_names = T
	if ( is.null(full_report) ) { 
		full_report <- a_report[,c(5,4),drop=F]
		report_md <- a_report[,c(5,3,2,1)]
	} else {
		full_report <- merge(full_report, a_report[,c(5,4),drop=F], by.x=1, by.y=1, all=T)
		report_md <- merge(report_md, a_report[,c(5,3,2,1)], by.x=1, by.y=1, all=T)
		report_md[is.na(report_md[,2]),2] <- report_md[is.na(report_md[,2]),4]
		report_md[is.na(report_md[,3]),3] <- report_md[is.na(report_md[,3]),5]
		report_md[is.na(report_md[,4]),4] <- report_md[is.na(report_md[,4]),6]
		report_md <- report_md[,1:4]
	}

	# Create a combined report of all samples
	last_col <- ncol(ecm_data)
	colnames(ecm_data)[7:last_col] <- paste(sample_name, colnames(ecm_data)[7:last_col], sep=" : ")
	ecm_data[,'Division'] <- as.character(ecm_data[, 'Division'])
	ecm_data[,'Category'] <- as.character(ecm_data[, 'Category'])
	if ( is.null(sample_details) ) { 
		sample_details <- ecm_data	
	} else {
		sample_details <- merge(sample_details, ecm_data[,c(1,7:last_col)], by.x=1, by.y=1, all=T)
		# Check if any proteins are missing annotations...
		row.names(ecm_data) <- as.character(ecm_data$ProteinID)
		row.names(sample_details) <- as.character(sample_details$ProteinID)
		missing_ids <- row.names(sample_details)[is.na(sample_details$Division)]
		# For any IDs with missing annotations, add from the current ecm_data object
		for ( an_id in missing_ids ) { 
			sample_details[an_id, c('Division', 'Category', "Gene symbol", 'Gene name', 'Species')] <- 
				ecm_data[an_id, c('Division', 'Category', "Gene symbol", 'Gene name', 'Species')]
		}
	}
}



# If specified, write out the combined protein report
if ( ! is.null(opts$details) ) { 
	sample_details[is.na(sample_details)] <- 0
	write.table(sample_details[order(sample_details[,1]), ], sep="\t", file=opts$details, row.names=F, quote=F)
}

ecm_data[,'Division'] <- factor(ecm_data[, 'Division'], levels=c("Core matrisome", "Matrisome-associated", "Other"))
ecm_data[,'Category'] <- factor(ecm_data[, 'Category'])

colnames(report_md) <- sub(".x", "", colnames(report_md), fixed=T)

row.names(full_report) <- full_report$column
# Delete the first column, helps with cbinding
full_report <- full_report[,-1, drop=F]
# Set any NA values to zero
full_report[is.na(full_report)] <- 0

full_report <- full_report[ ! grepl("^Percent.+:TOTAL$", row.names(full_report)), , drop=F]

full_report <- t(full_report)
full_report <- data.frame(Sample=row.names(full_report), full_report, check.names=F)

# Write the report to a file
if ( grepl("\\.xlsx$", opts$output) ) { 
	# If the report file is an xlsx file, then generate a Excel file with tabs for each section
	library(xlsx)
	wb <- createWorkbook()

	# Setup variables to format the columns with percent values
	percent_format <- CellStyle(wb, dataFormat=DataFormat("0.00%"))
	percent_start <- length(grep("^Protein_count", colnames(full_report))) + 2
	percent_end <- length(grep("^Percent_protein", colnames(full_report))) + percent_start
	percent_cols <- percent_start:percent_end
	table_styles <- list()
	for ( i in 1:length(percent_cols) ) { 
		table_styles[[i]] <- percent_format
	}
	names(table_styles) <- as.character(percent_cols)

	# Create sheet for protein counts
	sheet <- createSheet(wb, sheetName="Protein count")
	addDataFrame(full_report[,c(1, grep("^Protein_count", colnames(full_report)), grep("^Percent_protein", colnames(full_report)))], sheet,
		row.names=F, colStyle=table_styles)

	# Create sheet for peptide counts
	sheet <- createSheet(wb, sheetName="Peptide count")
	addDataFrame(full_report[,c(1, grep("^Peptide_count", colnames(full_report)), grep("^Percent_peptides", colnames(full_report)))], sheet,
		row.names=F, colStyle=table_styles)

	# Create sheet for spectra counts
	sheet <- createSheet(wb, sheetName="Spectra count")
	addDataFrame(full_report[,c(1, grep("^Spectra_count", colnames(full_report)), grep("^Percent_spectra", colnames(full_report)))], sheet,
		row.names=F, colStyle=table_styles)

	# Create sheet for EIC intensity
	sheet <- createSheet(wb, sheetName="EIC intensity")
	addDataFrame(full_report[,c(1, grep("^EIC_intensity", colnames(full_report)), grep("^Percent_EIC", colnames(full_report)))], sheet,
		row.names=F, colStyle=table_styles)

	# Save the workbook to an Excel file
	saveWorkbook(wb, file=opts$output)
} else {
	# Otherwise create a tab delimited report
	# First re-arrange the columns to be protein, peptide, spectra, and EIC	
	full_report <- full_report[,c(1, grep("^Protein_count", colnames(full_report)), grep("^Percent_protein", colnames(full_report)),
		grep("^Peptide_count", colnames(full_report)), grep("^Percent_peptides", colnames(full_report)),
		grep("^Spectra_count", colnames(full_report)), grep("^Percent_spectra", colnames(full_report)),
		grep("^EIC_intensity", colnames(full_report)), grep("^Percent_EIC", colnames(full_report)))]
	write.table(full_report, sep="\t", file=opts$output, row.names=F, quote=F)
}

##
# function to quickly prep the plot data from the full report
prepPlotData <- function(x, value.name, divisions=c()) { 
	max_col <- ncol(x)
	# sel_cols <- c("Sample", grep(paste0("^", value.name), colnames(x), value=T))
	plot_data <- x[ x$variable == value.name, 3:max_col, drop=F]
	plot_data <- melt(plot_data, id.vars=c("Species","Group"), variable.name="Sample", value.name=value.name)		
	plot_data <- plot_data[ plot_data$Group %in% divisions, ]	
	plot_data[,value.name] <- as.numeric(plot_data[,value.name])
	return(plot_data)
}

### Generate the plots
if ( ! is.null(opts$plot) ) { 
	library(ggplot2)
	# Setup the plot
	png_plot = ( opts$type == "png" )
	if ( ! png_plot ) { pdf(opts$plot, height=opts$height, width=opts$width) } 
	divisions <- levels(ecm_data$Division)
	# Remove the total columns
	report_data <- merge(report_md, t(full_report), by.x=1, by.y=0)
	report_data <- report_data[ report_data$Group != "TOTAL", ]
	# base style for pie charts
	if ( opts$style == "pie" ) {
		pie_style <- theme_minimal() + 
			theme(axis.title.x = element_blank(), 
				axis.title.y = element_blank(), 
				panel.border = element_blank(), 
				panel.grid=element_blank(), 
				axis.text.x=element_blank(),
				axis.text.y=element_blank(),
				axis.ticks = element_blank(),
				text = element_text(size=16))
	}

	facet_cols = ifelse(length(unique(counts_data$Sample)) > 10, 
		ifelse(length(unique(counts_data$Sample)) > 50, 10, 5), 3)

	# Histograms of counts
	if ( png_plot ) { png(paste0(opts$plot, ".hist_proteins.png"), width=opts$width, height=opts$height) }
	print(ggplot(counts_data, aes(Peptides, fill=Division)) + geom_histogram(bins=30) + xlab("Peptides") + 
		ylab("Count of proteins") + ggtitle("Distribution of peptide counts") + scale_fill_manual(values = matrisome_palette) +
		theme_bw() + facet_wrap(. ~ Sample, ncol=facet_cols))

	if ( png_plot ) { 
		invisible(dev.off())
		png(paste0(opts$plot, ".hist_psms.png"), width=opts$width, height=opts$height) 
	}
	print(ggplot(counts_data, aes(Spectra, fill=Division)) + geom_histogram(bins=30) + xlab("PSMs") + 
		ylab("Count of proteins") + ggtitle("Distribution of Spectra/PSM counts") + scale_fill_manual(values = matrisome_palette) +
		theme_bw() + facet_wrap(. ~ Sample, ncol=facet_cols))

	# Protein plot
	if ( png_plot ) { 
		invisible(dev.off())
		png(paste0(opts$plot, ".proteins.png"), width=opts$width, height=opts$height) 
	}

	species_list <- levels(report_data$Species)
	species_list <- species_list[ species_list != "ALL" ]

	has_species <- length(species_list) > 1

	if ( ! has_species ) { 
		report_data <- report_data[ report_data$Species != "ALL", ]
	}

	plot_data <- prepPlotData(report_data, "Protein_count", divisions)
	colnames(plot_data)[2] <- "Division"
	if ( opts$style == "pie" ) { 
		baseplot <- ggplot(plot_data, aes(x=1, y=Protein_count, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle('Protein count') +
			scale_fill_manual(values = matrisome_palette) +
			coord_polar("y", start=0, direction=-1)
		if ( has_species ) { 
			print(baseplot + facet_grid(Sample ~ Species))
		} else {
			print(baseplot + facet_wrap(. ~ Sample, ncol=facet_cols))
		}
	} else {
		baseplot <- ggplot(plot_data, aes(x=Sample, y=Protein_count, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			geom_bar(stat='identity') + theme_bw()
		if ( has_species ) { 
			print(baseplot + facet_wrap(. ~ Species, nrow=1))
		} else {
			print(baseplot)
		}
	}

	# Peptide plot
	if ( png_plot ) { 
		invisible(dev.off())
		png(paste0(opts$plot, ".peptides.png"), width=opts$width, height=opts$height) 
	}
	plot_data <- prepPlotData(report_data, "Peptide_count", divisions)
	colnames(plot_data)[2] <- "Division"
	if ( opts$style == "pie" ) { 
		baseplot <- ggplot(plot_data, aes(x=1, y=Peptide_count, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("Peptide count") +
			scale_fill_manual(values = matrisome_palette) +
			coord_polar("y", start=0, direction=-1)
		if ( has_species ) { 
			print(baseplot + facet_grid(Sample ~ Species))
		} else {
			print(baseplot + facet_wrap(. ~ Sample, ncol=facet_cols))
		}
	} else {
		baseplot <- ggplot(plot_data, aes(x=Sample, y=Peptide_count, fill=Division)) + 
			geom_bar(stat='identity') + scale_fill_manual(values = matrisome_palette) +
			theme_bw() + ylab("peptide count")
		if ( has_species ) { 
			print(baseplot + facet_wrap(. ~ Species, nrow=1))
		} else {
			print(baseplot)
		}
	}

	# Spectra/PSM plot
	if ( png_plot ) { 
		invisible(dev.off())
		png(paste0(opts$plot, ".psms.png"), width=opts$width, height=opts$height) 
	}
	plot_data <- prepPlotData(report_data, "Spectra_count", divisions)
	colnames(plot_data)[2] <- "Division"
	if ( opts$style == "pie" ) { 
		baseplot <- ggplot(plot_data, aes(x=1, y=Spectra_count, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("Spectra count") +
			scale_fill_manual(values = matrisome_palette) +
			coord_polar("y", start=0, direction=-1)
		if ( has_species ) { 
			print(baseplot + facet_grid(Sample ~ Species))
		} else {
			print(baseplot + facet_wrap(. ~ Sample, ncol=facet_cols))
		}
	} else {
		baseplot <- ggplot(plot_data, aes(x=Sample, y=Spectra_count, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			geom_bar(stat='identity') + theme_bw()
		if ( has_species ) { 
			print(baseplot + facet_wrap(. ~ Species, nrow=1))
		} else {
			print(baseplot)
		}
	}

	# EIC plot
	if ( png_plot ) { 
		invisible(dev.off())
		png(paste0(opts$plot, ".eic.png"), width=opts$width, height=opts$height) 
	}
	plot_data <- prepPlotData(report_data, "EIC_intensity", divisions)
	colnames(plot_data)[2] <- "Division"
	if ( opts$style == "pie" ) { 
		baseplot <- ggplot(plot_data, aes(x=1, y=EIC_intensity, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("MS1 intensity") +
			scale_fill_manual(values = matrisome_palette) +
			coord_polar("y", start=0, direction=-1)
		if ( has_species ) { 
			print(baseplot + facet_grid(Sample ~ Species))
		} else {
			print(baseplot + facet_wrap(. ~ Sample, ncol=facet_cols))
		}
	} else {
		baseplot <- ggplot(plot_data, aes(x=Sample, y=EIC_intensity, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			facet_wrap(. ~ Species, nrow=1) +
			geom_bar(stat='identity') + theme_bw()
		if ( has_species ) { 
			print(baseplot + facet_wrap(. ~ Species, nrow=1))
		} else {
			print(baseplot)
		}
	}

	invisible(dev.off())
}
