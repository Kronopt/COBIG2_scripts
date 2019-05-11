#!/bin/bash

# Run Blast2GO for every FASTA file in a folder (recursively)
# Must be run inside the folder containing the blast2go_cli.run executable

# SCRIPT PARAMETERS
# $1 = Folder containing FASTA files
# $2 = Output folder
# $3 = Temporary blastXML folder

# All parameters must be supplied
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo "Must be run inside the folder containing the blast2go_cli.run executable"
	echo
        echo "3 arguments needed:"
        echo "1 - Folder containing FASTA files"
        echo "2 - Output folder"
        echo "3 - Temporary blastXML folder (must be defined in the blast2go properties file as well"
        exit 1
fi

# clean blastXML folder if it is not empty
if ! [ -z $(ls -A) ]; then
	rm $3*
fi

# run blast2go for each FASTA file
for fasta in $(find $1 -name '*.fasta'); do
	fasta_name=$(basename $fasta)
	echo
	echo "############"
	echo Blast2GO $fasta_name
	echo "############"
	echo

	./blast2go_cli.run -properties cli.prop -loadfasta $fasta -localblast binblast/ \
	-mapping -annotation -useobo go_latest.obo.gz \
	-savelog $2/$fasta_name/$fasta_name\_log.txt \
	-savereport $2/$fasta_name/$fasta_name\_report.pdf \
	-saveseqtable $2/$fasta_name/$fasta_name\_seqtable.tab \
	-statistics all -nameprefix $fasta_name -workspace $2

	# clean blastXML folder
	rm $3*
done

echo "done"
