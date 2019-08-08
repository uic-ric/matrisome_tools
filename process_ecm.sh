#!/bin/sh
################################################################################
# Script : ecm_process.sh
# Author : George Chlipala
# Created: Feb 13, 2019
# -- Description ----------------------------------------
# Master script to process mzML and pepXML file to obtain EIC intensities
# for each peptide then annotate proteins using ECM database
# -- Requirements ---------------------------------------
# OpenMS
################################################################################

if [ "$#" -lt "4" ] ; then
	echo "usage: $0 [mzML file] [idXML/pepXML/mzIdentML file] [ECM DB file] [output protein table] [sample name]"
	echo 
	echo "    Master script to process mzML and idXML/pepXML file to obtain PSM counts and EIC intensities"
	echo "    for each peptide then annotate proteins using ECM database"
	exit
fi

# Get the variables
mzfile="$1"
idxmlfile="$2"
ecm_db="$3"
# Output
protein_table="$4"

if [ "$#" -eq "5" ] ; then
	# Set the sample name from the variable
	sample_name="$5"
else
	# Set the sample name from the mzML file
	sample_name=$(basename $mzfile)
	sample_name=${sample_name%.*}
fi

temp_files=""

if [ "$ID_FILE_TYPE" != "" ] ; then
	# If running in galaxy, all files have a .dat extension.  
	# OpenMS requires that files have the proper extension, e.g. .pepXML,
	# So, create a symbolic link with the proper extension 
	case $ID_FILE_TYPE in 
		pepxml|protXML|mascotXML|omssaXML|idXML)
			ID_FILE_TYPE=$(echo $ID_FILE_TYPE | sed -e 's/xml/XML/') ;;
	esac
	in_file="$(mktemp -u tmp_ecm_process_convertXXXXXX).$ID_FILE_TYPE"
	ln -s $idxmlfile $in_file
	temp_files="$in_file $temp_files"
	idxmlfile="$in_file"
fi

if [ "$MS_FILE_TYPE" != "" ] ; then
	# If running in galaxy, all files have a .dat extension.  
	# OpenMS requires that files have the proper extension, e.g. .pepXML,
	# So, create a symbolic link with the proper extension 
	case $MS_FILE_TYPE in 
		mzml) MS_FILE_TYPE="mzML" ;;
		mzxml) MS_FILE_TYPE="mzXML" ;;
		mzdata) MS_FILE_TYPE="mzData" ;;
	esac
	new_file="$(mktemp -u tmp_ecm_process_convertXXXXXX).$MS_FILE_TYPE"
	ln -s $mzfile $new_file
	mzfile="$new_file"
	temp_files="$new_file $temp_files"
fi

# check if the ID file is an idXML file.  If not, convert using IDFileConverter
if [ "${idxmlfile##.*}" != "idXML" ] ; then
	converted_file="$(mktemp -u tmp_ecm_process_convertXXXXXX).idXML"
	echo Running IDFileConverter 
	if [ "${idxmlfile##.*}" != "pepXML" ] ; then
		mz_name=$(grep '<msms_run_summary' ${idxmlfile} | sed -e 's/^.*base_name="//' | cut -d\" -f1)
		#mz_name=$(grep search_summary ${idxmlfile} | sed -e 's/^.+base_name="(.+)".+$/\\1/')
		if [ "$mz_name" == "" ] ; then
			pepxmlfile="$(mktemp -u tmp_ecm_process_convertXXXXXX).pepXML"
			sed -e 's/ base_name=""/ base_name="ECM"/' < $idxmlfile > $pepxmlfile
			mz_name="ECM"
			idxmlfile="$pepxmlfile"
			temp_files="$pepxmlfile $temp_files"
		fi
		conv_args="-mz_name $mz_name"
	fi
	echo IDFileConverter -in $idxmlfile -out $converted_file -mz_file $mzfile -no_progress $conv_args
	IDFileConverter -in $idxmlfile -out $converted_file -mz_file $mzfile -no_progress $conv_args 2>&1
	if [ $? != 0 ] ; then
		rc=$?
		>&2 echo "ERROR: Could not convert input ID file"
		exit $rc
	fi
	idxmlfile=$converted_file
fi

## WORKAROUND
# Filter file to remove proteins with ambiguous amino acids
filtered_file="$(mktemp -u tmp_ecm_process_filteredXXXXXX).idXML"
echo IDFilter -in $idxmlfile -out $filtered_file -blacklist:protein_accessions Q7M0E4 A0A088QCC0
IDFilter -in $idxmlfile -out $filtered_file -blacklist:protein_accessions Q7M0E4 A0A088QCC0 2>&1
if [ $? != 0 ] ; then
	rc=$?
	>&2 echo "ERROR: Could not filter input ID file"
	exit $rc
fi
idxmlfile=$filtered_file

# Generate feature file to contain EIC intensities for all peptides
feature_xml="$(mktemp -u tmp-ecm_process-featureXXXXXX).featureXML"
temp_files="$feature_xml $temp_files"
echo FeatureFinderIdentification -in $mzfile -id $idxmlfile -out $feature_xml
FeatureFinderIdentification -in $mzfile -id $idxmlfile -out $feature_xml -no_progress 2>&1
if [ $? != 0 ] ; then
	rc=$?
	>&2 echo "ERROR: Could detect protein features"
	exit $rc
fi

# Generate tables for PSM counts
psm_protein="$(mktemp -u tmp-ecm_process-psm_proteinXXXXXX).csv"
temp_files="$psm_protein $temp_files"
ProteinQuantifier -in $idxmlfile -out $psm_protein -include_all -average sum 2>&1
if [ $? != 0 ] ; then
	rc=$?
	>&2 echo "ERROR: Could quantitate PSMs"
	exit $rc
fi

# Generate tables for EIC counts
eic_protein="$(mktemp -u tmp-ecm_process-eic_proteinXXXXXX).csv"
temp_files="$eic_protein $temp_files"
ProteinQuantifier -in $feature_xml -out $eic_protein -include_all -average sum
if [ $? != 0 ] ; then
	rc=$?
	>&2 echo "ERROR: Could quantitate peptide EICs"
	exit $rc
fi

# Run script to merge the results
python ${ECM_PY_SCRIPT:=process_ecm.py} -p $psm_protein -e $eic_protein --db $ecm_db -o $protein_table -s "$sample_name"

# Cleanup temporary files
for tmp in "$tempfiles" ; do
	rm -f $tmp
done
