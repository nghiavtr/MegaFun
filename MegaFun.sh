#!bin/bash


syntax="/path/to/MegaFun/MegaFun.sh \n\t-in [input folder with gz files] \n\t-pan [pangenome in fasta format] \n\t-iso [isolate genomes in fasta format] \n\t-smap [species mapping between chocophlan and ncbi] \n\t-p [number of thread] "

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

    -in|--input) 
    sample_Dir=$2
    shift
    ;; 

    -pan|--pangenome) 
    pan_Dir=$2
    shift
    ;; 


    -iso|--isolate) 
    iso_Dir=$2
    shift
    ;; 

    -smap|--speciesmap) 
    smapFile=$2
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


if [[ -z "$sample_Dir" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input folder provided !"
   echo ""
   exit
fi

if [[ -z "$pan_Dir" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no pan genome input folder !"
   echo ""
   exit
fi


if [[ -z "$iso_Dir" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no isolate genome input folder !"
   echo ""
   exit
fi

if [[ -z "$smapFile" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no mapping specices information (between chocophlan and ncbi), please run extract_isolate_list.R to obtain it !"
   echo ""
   exit
fi



if [[ -z "$ncore" ]]; then
   ncore=8
   echo "No specific setting for thread number, so use the default setting (-p 8)"
fi




pan_eqv_Dir="pan_eqv"

#clean up temp folders
if [[ -d temp ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ temp/
    rm -rf empty_dir00
    rm -rf temp
fi
if [[ -d sam_temp ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ sam_temp/
    rm -rf empty_dir00
    rm -rf sam_temp
fi

if [[ -d $pan_eqv_Dir ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ $pan_eqv_Dir
    rm -rf empty_dir00
    rm -rf $pan_eqv_Dir
fi




genomeFastaFile="pangenome_combined.fa"
cat $pan_Dir/*.ffn > $genomeFastaFile


bowtie2_indexFolder="bowtie2_ref"
bowtie2_index=$(echo $bowtie2_indexFolder"/genome")



###create bowtie2_ref, the index folder of bowtie2 for genome.fa (the pangenome fasta sequence of species) 
echo "Create bowtie2 index"
mkdir $bowtie2_indexFolder
bowtie2-build --threads $ncore $genomeFastaFile $bowtie2_index

#### design X matrix #####
echo "Generate X_matrix"

# simulate for X
mkdir X_matrix
/path/to/MegaFun/bin/simx $pan_Dir -o X_matrix/temp -t $ncore
cat X_matrix/temp/*.* > X_matrix/sample_01.fasta.gz
rm X_matrix/temp -f -r

# mapping for X
#minimap2 --secondary=yes -t $ncore -Y -ax sr $genomeFastaFile X_matrix/sample_01.fasta.gz | samtools sort -n -o X_matrix/sample_01.bam
minimap2 --secondary=yes -t $ncore -Y -ax sr $genomeFastaFile X_matrix/sample_01.fasta.gz | samtools view -b -o  X_matrix/sample_01.bam


# generate eqClasses for X
/path/to/MegaFun/bin/gentc X_matrix/sample_01.bam -o X_matrix/eqClass.txt -nm 5 -t $ncore

# create X
Rscript /path/to/MegaFun/buildCRP.R


# clean simulated data 
rm X_matrix/sample_01.bam
rm X_matrix/sample_01.fasta.gz

###### generate Y count #####
echo "Create Y count"
#create output folder for eqclasses
mkdir $pan_eqv_Dir

for FILE in `ls ${sample_Dir}/*.gz`
do
    echo "${FILE}"
    prefix=`basename ${FILE} | cut -d. -f1`
    in_fq=$(echo "${FILE}")
    out_bam=$(echo "${pan_eqv_Dir}/${prefix}.bam")
    out_bam_sorted=$(echo "${pan_eqv_Dir}/${prefix}_sort.bam")
    out_eqc_folder=$(echo "${pan_eqv_Dir}/${prefix}")
    out_eqc=$(echo $out_eqc_folder"/eqClass.txt")
    echo "$in_fq"
    echo "$out_bam"
    echo "$out_bam_sorted"
    echo "$out_eqc"
    echo -e "\n"

    # create a folder for the sample
    mkdir $out_eqc_folder

    # mapping for data
    bowtie2 -f -x $bowtie2_index -U $in_fq -p $ncore --very-sensitive -a | samtools sort -n -o $out_bam

    # eqClasses for data 
    /path/to/MegaFun/bin/bam2eqv $out_bam -o $out_eqc -t 1 -nm_nonsharing 30
    #sort the bam file
    samtools sort $out_bam -o $out_bam_sorted
    samtools index $out_bam_sorted
    rm $out_bam
done



# generate Y count
Rscript /path/to/MegaFun/create_count_matrix.R in=${pan_eqv_Dir} > create_count_matrix.log 2>&1

##### AEM for Xmatrix from pan-genome #####
echo "AEM for pan-X"
Rscript /path/to/MegaFun/AEM_update_X_beta.R in=${pan_eqv_Dir} > AEM_update_X_beta.log 2>&1 

# Re-design X matrix with isolate reference
echo "Redesign for isolate-X"

# step 1: create X matrix - it is independent from step 2 - work on temp folder
#create the map between the short name and original name of pangenome
Rscript /path/to/MegaFun/panNameMap.R pangenome=$genomeFastaFile > panNameMap.log 2>&1
echo "split isolate fasta sequences into CRPs"
Rscript /path/to/MegaFun/simx_isolate_splitIsoSeq.R iso=$iso_Dir pangenome=$genomeFastaFile smap=$smapFile >out_iso_splitIsoSeq.log 2>&1
Rscript /path/to/MegaFun/simx_isolate_addPanGenSeq.R pangenome=$genomeFastaFile >out_iso_addPanGenSeq.log 2>&1

echo "build X matrix at isolate level"
Rscript /path/to/MegaFun/simx_isolate_buildX.R >out_iso_buildX.log 2>&1

## step 2: split bam files into corresponding CRPs. This will scan each bam file only once. It is independent from step 1
# perform all bam files in parallel at once
echo "split bam files into CRPs"
mkdir sam_temp
Rscript /path/to/MegaFun/simx_isolate_splitBam.R in=${pan_eqv_Dir} pangenome=$genomeFastaFile >out_iso_splitBam.log 2>&1

## step 3: create Y matrix - it requires outputs of step 1 and step 2
mkdir output
echo "create Y count at isolate level"
Rscript /path/to/MegaFun/simx_isolate_createY.R in=${pan_eqv_Dir} >out_iso_createY.log 2>&1
# Need to reuse some Y from pan-genome due to no exist isolate sequences
Rscript /path/to/MegaFun/combine_Y.R pangenome=$genomeFastaFile in=${pan_eqv_Dir} > combine_Y.log 2>&1 

# Rerun AEM for new isolate-X matrix
echo "AEM for isolate-X"
Rscript /path/to/MegaFun/AEM_isolate.R > AEM_isolate.log 2>&1 

# Copy quantification results at pangenome level to the ouput
mkdir output/Pan
cp ${pan_eqv_Dir}/*.RData output/Pan

## step 4: clean up intermediate data - activate this after we succesfully tested the tool to release storage space
if [[ -d temp ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ temp/
    rm -rf empty_dir00
    rm -rf temp
fi
if [[ -d sam_temp ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ sam_temp/
    rm -rf empty_dir00
    rm -rf sam_temp
fi
if [[ -d $pan_eqv_Dir ]]; then
    mkdir empty_dir00
    rsync -a --delete empty_dir00/ $pan_eqv_Dir
    rm -rf empty_dir00
    rm -rf $pan_eqv_Dir
fi


# output:
## 1: final_paralog.RData
## 2: final_singleton.RData
