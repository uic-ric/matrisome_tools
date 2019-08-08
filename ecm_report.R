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
library(ggplot2)

##
# Setup script arguments
option_list <- list( 
    make_option(c("-o", "--output"), help="Output report (text). Default is STDOUT.", default="-"),
    make_option(c("-d", "--details"), help="Output details for each sample (text)."),
    make_option(c("-m", '--min_peptides'), help="Minimum number of peptides. Default is %default", type="integer", default=2),
    make_option(c("-p", "--plot"), help="Output plots."),
    make_option(c("-t", "--type"), help="Output type, either png or pdf. Default is %default.", default="pdf"),
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
	positional_arguments = c(1,Inf) )

# Get options
opts <- args$options

# If help was specified, print usage and quit
if ( opts$help ) { printHelp(opt_parser) }

# Parse input files
data_files = character()

for ( arg in args$args ) { 
	files = Sys.glob(arg)
	if ( length(files) > 0 ) { 
		data_files = c(data_files, files)
	} else {
		data_files = c(data_files, arg)
	}
}

##
# function to parse each sample file
# Arguments are..
# 1. data frame, of the input data
# 2. list of categories to include. Default is Collagens
### Columns of input data 
### Protein_ID	Division	Category	Peptides_spectra	Peptides_EIC	PSM_count	EIC_intensity 
sampleReport <- function(x, categories=c("Collagens")) { 
	# Generate data frame to contain report
	report <- as.data.frame( matrix(nrow=(length(levels(x$Division)) + length(categories) + 1), ncol=9) )
	row.names(report) <- c( levels(x$Division), categories, "TOTAL")
	colnames(report) <- c("Protein_count", "Percent_protein", "Peptides_spectra", "Percent_spectra", "PSM_count", "Percent_PSM", "Peptides_EIC", "EIC_intensity", "Percent_EIC")

	# First, generate the totals and relative stats
	report["TOTAL", "Protein_count"] <- nrow(x)
	if ( "Peptides_spectra" %in% colnames(x) ) { 
		report["TOTAL", "Peptides_spectra"] = sum(x$Peptides_spectra)
		report["TOTAL", "PSM_count"] = sum(x$PSM_count)
	}
	if ( "Peptides_EIC" %in% colnames(x) ) { 
		report["TOTAL", "Peptides_EIC"] = sum(x$Peptides_EIC)
		report["TOTAL", "EIC_intensity"] = sum(x$EIC_intensity)
	}

	# Next, generate the report for the major divisions
	for ( division in levels(x$Division) ) { 
		if ( "Peptides_spectra" %in% colnames(x) ) { 
			report[division, "Peptides_spectra"] <- sum(x[x$Division == division, "Peptides_spectra"])
			report[division, "PSM_count"] <- sum(x[x$Division == division, "PSM_count"])
			report[division, "Percent_spectra"] <- report[division, "Peptides_spectra"] / report["TOTAL", "Peptides_spectra"]
			report[division, "Percent_PSM"] <- report[division, "PSM_count"] / report["TOTAL", "PSM_count"]
		}
		if ( "Peptides_EIC" %in% colnames(x) ) { 
			report[division, "Peptides_EIC"] <- sum(x[x$Division == division, "Peptides_EIC"])
			report[division, "EIC_intensity"] <- sum(x[x$Division == division, "EIC_intensity"])
			report[division, "Percent_EIC"] <- report[division, "EIC_intensity"] / report["TOTAL", "EIC_intensity"]
		}
		report[division,"Protein_count"] <- length(which(x$Division == division))
		report[division, "Percent_protein"] <- report[division, "Protein_count"] / report["TOTAL", "Protein_count"]
	}	

	# Finally, generate the report for the selected categories
	for ( c in categories ) { 
		if ( "Peptides_spectra" %in% colnames(x) ) { 
			report[c, "Peptides_spectra"] <- sum(x[x$Category == c, "Peptides_spectra"])
			report[c, "PSM_count"] <- sum(x[x$Category == c, "PSM_count"])
			report[c, "Percent_spectra"] <- report[c, "Peptides_spectra"] / report["TOTAL", "Peptides_spectra"]
			report[c, "Percent_PSM"] <- report[c, "PSM_count"] / report["TOTAL", "PSM_count"]
		}
		if ( "Peptides_EIC" %in% colnames(x) ) { 
			report[c, "Peptides_EIC"] <- sum(x[x$Category == c, "Peptides_EIC"])
			report[c, "EIC_intensity"] <- sum(x[x$Category == c, "EIC_intensity"])
			report[c, "Percent_EIC"] <- report[c, "EIC_intensity"] / report["TOTAL", "EIC_intensity"]
		}
		report[c,"Protein_count"] <- length(which(x$Category == c))
		report[c, "Percent_protein"] <- report[c, "Protein_count"] / report["TOTAL", "Protein_count"]
	}
	
	# Return the report
	return(report)
}

# Parse over the different samples and generate the full report
full_report = NULL
sample_details = NULL
counts_data = data.frame(Sample=character(), Division=character(), Protein_ID=character(), Peptides_spectra=numeric(), PSM_count=numeric())

for ( sample_file in data_files ) { 
	## Read inputfile
	ecm_data <- read.delim(sample_file, comment.char='#', check.names=F)
	ecm_data[,'Gene name'] <- as.character(ecm_data[, 'Gene name'])
	## Read the comment lines, if they exist
	ecm_comments <- grep("^##", readLines(sample_file, n=10), value=T)
	if ( length(ecm_comments) > 0 ) { 
		for ( line in ecm_comments ) { 	
			parts <- unlist(strsplit(line, ':'))		
			field <- tolower(trimws(sub("^##+ *", "", parts[1])))
			value <- trimws(paste(parts[-1], collapse=':'))
			if ( field == "sample" ) { sample_name = value }
		}
	} else {
		sample_name = sub("\\.(txt|tsv|csv|tab)$", "", basename(sample_file))
	}

	ecm_data <- ecm_data[ecm_data$Peptides_spectra >= opts$min_peptides, ]

	# Get counts for each protein.  For histograms in report
	sample_counts = data.frame(Sample=rep(sample_name, nrow(ecm_data)), ecm_data[,c("Protein_ID", "Division", "Peptides_spectra", "PSM_count")])
	counts_data = rbind(counts_data, sample_counts)
	
	# Generate the report for this sample
	a_report <- sampleReport(ecm_data)
	# Melt the report table
	a_report$columns <- row.names(a_report)
	a_report <- melt(a_report, id.vars="columns")
	row.names(a_report) <- paste(a_report$variable, a_report$columns, sep=":")		
	# Set the value column to the sample name
	colnames(a_report)[3] <- sample_name
	# Add this sample to the full report
	if ( is.null(full_report) ) { 
		full_report <- a_report[,3,drop=F]
	} else {
		full_report <- cbind(full_report, a_report[,3,drop=F])
	}

	# Create a combined report of all samples
	last_col <- ncol(ecm_data)
	colnames(ecm_data)[5:last_col] <- paste(sample_name, colnames(ecm_data)[5:last_col], sep=" : ")
	if ( is.null(sample_details) ) { 
		sample_details <- ecm_data	
	} else {
		sample_details <- merge(sample_details, ecm_data[,c(1,5:last_col)], by.x=1, by.y=1, all=T)
		# Check if any proteins are missing annotations...
		row.names(ecm_data) <- ecm_data$Protein_ID
		row.names(sample_details) <- sample_details$Protein_ID
		missing_ids <- row.names(sample_details)[is.na(sample_details$Division)]
		# For any IDs with missing annotations, add from the current ecm_data object
		for ( an_id in missing_ids ) { 
			sample_details[an_id, c('Division', 'Category', 'Gene name')] <- 
				ecm_data[an_id, c('Division', 'Category', 'Gene name')]
		}
	}
}

# If specified, write out the combined protein report
if ( ! is.null(opts$details) ) { 
	write.table(sample_details[order(sample_details[,1]), ], sep="\t", file=opts$details, row.names=F, quote=F)
}

full_report <- full_report[ ! grepl("^Percent.+:TOTAL$", row.names(full_report)), , drop=F]

full_report <- t(full_report)
full_report <- data.frame(Sample=row.names(full_report), full_report, check.names=F)

# Write the report to a file
write.table(full_report, sep="\t", file=opts$output, row.names=F, quote=F)

##
# function to quickly prep the plot data from the full report
prepPlotData <- function(x, value.name, divisions=c()) { 
	sel_cols <- c("Sample", grep(paste0("^", value.name), colnames(x), value=T))
	plot_data <- x[, sel_cols, drop=F]
	plot_data <- melt(plot_data, id.vars="Sample", variable.name="Division", value.name=value.name)		
	plot_data$Division <- sub("^.+:", "", plot_data$Division)
	plot_data <- plot_data[ plot_data$Division %in% divisions, ]	
	return(plot_data)
}

### Generate the plots
if ( ! is.null(opts$plot) ) { 
	pdf(opts$plot)
	divisions <- levels(ecm_data$Division)
	# Remove the total columns
	report_data <- full_report[, ! grepl("TOTAL$", colnames(full_report)), drop=F]	
	# base style for pie charts
	if ( opts$style == "pie" ) {
		pie_style <- theme_minimal() + 
			theme(axis.title.x = element_blank(), 
				axis.title.y = element_blank(), 
				panel.border = element_blank(), 
				panel.grid=element_blank(), 
				axis.text.x=element_blank(),
				axis.text.y=element_blank(),
				axis.ticks = element_blank())
	}

	# Histograms of counts
	print(ggplot(counts_data, aes(Peptides_spectra, fill=Division)) + geom_histogram(bins=30) + xlab("Peptides") + 
		ylab("Count of proteins") + ggtitle("Distribution of peptide counts") + scale_fill_manual(values = matrisome_palette) +
		theme_bw() + facet_wrap(. ~ Sample, ncol=3))

	print(ggplot(counts_data, aes(PSM_count, fill=Division)) + geom_histogram(bins=30) + xlab("PSMs") + 
		ylab("Count of proteins") + ggtitle("Distribution of PSM counts") + scale_fill_manual(values = matrisome_palette) +
		theme_bw() + facet_wrap(. ~ Sample, ncol=3))

	# Protein plot
	plot_data <- prepPlotData(report_data, "Protein_count", divisions)
	if ( opts$style == "pie" ) { 
		print(ggplot(plot_data, aes(x=1, y=Protein_count, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle('Protein count') +
			scale_fill_manual(values = matrisome_palette) +
			facet_wrap(. ~ Sample, ncol=3) + coord_polar("y", start=0, direction=-1))
	} else {
		print(ggplot(plot_data, aes(x=Sample, y=Protein_count, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			geom_bar(stat='identity') + theme_bw())
	}

	# Peptide plot
	plot_data <- prepPlotData(report_data, "Peptides_spectra", divisions)
	if ( opts$style == "pie" ) { 
		print(ggplot(plot_data, aes(x=1, y=Peptides_spectra, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("Peptide count") +
			scale_fill_manual(values = matrisome_palette) +
			facet_wrap(. ~ Sample, ncol=3) + coord_polar("y", start=0, direction=-1))
	} else {
		print(ggplot(plot_data, aes(x=Sample, y=Peptides_spectra, fill=Division)) + 
			geom_bar(stat='identity') + scale_fill_manual(values = matrisome_palette) +
			theme_bw() + ylab("peptide count"))
	}

	# Peptide plot
	plot_data <- prepPlotData(report_data, "PSM_count", divisions)
	if ( opts$style == "pie" ) { 
		print(ggplot(plot_data, aes(x=1, y=PSM_count, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("Spectra count") +
			scale_fill_manual(values = matrisome_palette) +
			facet_wrap(. ~ Sample, ncol=3) + coord_polar("y", start=0, direction=-1))
	} else {
		print(ggplot(plot_data, aes(x=Sample, y=PSM_count, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			geom_bar(stat='identity') + theme_bw())
	}

	# EIC plot
	plot_data <- prepPlotData(report_data, "EIC_intensity", divisions)
	if ( opts$style == "pie" ) { 
		print(ggplot(plot_data, aes(x=1, y=EIC_intensity, fill=Division)) + 
			geom_bar(stat='identity', position="fill") + pie_style + ggtitle("MS1 intensity") +
			scale_fill_manual(values = matrisome_palette) +
			facet_wrap(. ~ Sample, ncol=3) + coord_polar("y", start=0, direction=-1))
	} else {
		print(ggplot(plot_data, aes(x=Sample, y=EIC_intensity, fill=Division)) + 
			scale_fill_manual(values = matrisome_palette) +
			geom_bar(stat='identity') + theme_bw())
	}

	invisible(dev.off())
}
