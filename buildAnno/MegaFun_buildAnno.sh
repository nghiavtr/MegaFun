#!bin/bash

syntax="/path/to/MegaFun/buildAnno/MegaFun_buildAnno.sh \n\t-in [list of species] \n\t-out [output folder]  \n\t-pan [path to the folder of chocophlan pangenome] \n\t-map [mapping file to map between chocophlan and ncbi] \n\t-taxname [ncbi name mapping db] \n\t-taxdb [gz file of the summary of ncbi] \n\t-thread [number of thread] "

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

    -in|--input) 
    speciesList=$2
    shift
    ;;

    -out|--ouput) 
    anno_outpath=$2
    shift
    ;;  

    -pan|--pangenome) 
    chocoPath=$2
    shift
    ;;

    -map|--mapchoco) 
    chocoMap=$2
    shift
    ;;

    -taxname|--taxncbinamedb) 
    ncbi_name=$2
    shift
    ;; 
    

    -taxdb|--taxncbidb)    
    ncbiSummary=$2
    shift
    ;; 

    -p|--thread) 
    ncore=$2
    shift
    ;; 

    *)

esac
shift
done


if [[ -z "$speciesList" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no species list provided !"
   echo ""
   exit
fi


if [[ -z "$anno_outpath" ]]; then
   anno_outpath=$PWD
   echo "No specific output folder provided, the output will be save to the current folder!"
fi


if [[ -z "$chocoPath" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no pan genome input folder !"
   echo ""
   exit
fi


if [[ -z "$chocoMap" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no mapping file of pan genome !"
   echo ""
   exit
fi

if [[ -z "$ncbi_name" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no ncbi_name file !"
   echo ""
   exit
fi




if [[ -z "$ncbiSummary" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no ncbi summary file !"
   echo ""
   exit
fi


if [[ -z "$ncore" ]]; then
   ncore=8
   echo "No specific setting for thread number, so use the default setting (-p 8)"
fi


###################

#get the list of isolates
echo "Extract the list of isolates"
Rscript /path/to/MegaFun/buildAnno/extract_isolate_list.R speciesList=$speciesList out=$anno_outpath ncbiSummary=$ncbiSummary ncbiName=$ncbi_name chocoMap=$chocoMap

#download isolate annotation and extract pangenome sequence
echo "Download isolate annotation and extract pangenome sequence"
Rscript /path/to/MegaFun/buildAnno/get_pan_isolate_anno.R in=$anno_outpath/isolate_genome_summary.txt out=$anno_outpath/data chocoPath=$chocoPath 

#check if the download files are correct - redownload or clean
echo "Check and re-download"
Rscript /path/to/MegaFun/buildAnno/check_isolate_anno.R in=$anno_outpath/isolate_genome_summary.txt out=$anno_outpath/data chocoPath=$chocoPath 

# get CDS of genes in isolate, produce isolate_genome_summary.txt.filtered.csv to flag missed isolates
echo "Get CDS of isolate"
Rscript /path/to/MegaFun/buildAnno/get_isolate_CDS.R in=$anno_outpath/isolate_genome_summary.txt out=$anno_outpath/data

#vsearch to align isolate to pan to get protein id: This step is very slow due to vsearch performance. Consider to use map_isolate_proteins_chunk.R instead, where the number of rows in isolate_genome_summary.txt.filtered.csv is divided into N chunks (batches). Then we can submit each chunk separately for improving parallel process, alternatively you can use jobarray as instructed here https://www.pdc.kth.se/support/documents/run_jobs/job_arrays.html
echo "Align isolate to pan to get protein id"
Rscript /path/to/MegaFun/buildAnno/map_isolate_proteins.R in=$anno_outpath/isolate_genome_summary.txt.filtered.csv out=$anno_outpath/data

#get isolate genome
echo "Extract aligned gene sequence from isolate for MegaFun"
Rscript /path/to/MegaFun/buildAnno/get_isolate_genome.R in=$anno_outpath/isolate_genome_summary.txt.filtered.csv out=$anno_outpath/data

#shorten gene names in pangenome of chocophlanto avoid errors in samtools
Rscript /path/to/MegaFun/buildAnno/cut_genename_pan.R panIn=$anno_outpath/data/pangenome panOut=$anno_outpath/data/pangenome_shortNames

## clean up
echo "Clean up temp data"
if [[ -d $anno_outpath/data/isolate_genome ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ $anno_outpath/data/isolate_genome
    rm -rf empty_dir00
    rm -rf $anno_outpath/data/isolate_genome
fi