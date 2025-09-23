# MegaFun

MegaFun is the method for quantification of gene abundance from metagenomics data. Different from most existing tools, MegaFun considers i) multi-map reads between genes and ii) isolate-level sequences of species to improve the accuracy of the quantification.

## Prerequisite

- minimap2
- bowtie2
- samtools
- vsearch
- R and following R-packages: data.table, stringr, Biostrings, BSgenome, Rsamtools, parallel,foreach, doParallel

## Installation
- Download MegaFun, then configure it
```
bash config.sh
```
## How to run
### If the annotation is not prebuilt, need to create MegaFun annotation:

```
bash /path/to/MegaFun/buildAnno/MegaFun_anno_pipeline.sh -param params.txt
```

The output inludes species_list_summary.txt and two folders: isolate_genome_MegaFun and pangenome_shortNames, which will be used for the quantification step.

In the params.txt, following parameters are required

```
anno_outpath="/yourpath/annoOut"
chocoPath="/yourpath/full_chocophlan.v201901_v31/"
chocoMap="/yourpath/TaxChocoDB/choco_list_full_chocophlan.v201901_v31.RData"
ncbi_name="/yourpath/TaxChocoDB/ncbi_name_db.RData"
ncbiSummary="/yourpath/TaxChocoDB/assembly_summary_refseq.txt.gz"
ncore=16
metaphlan_outpath="/yourpath/metaphlan_out"
metaphlan_prop=0.5

runMetaPhlan="YES" # {YES, NO}
sample_Dir=/yourpath/Dataset # FASTA/FASTQ files for samples
suffix_sample="gz"
file_type="fastq" # {'fastq', 'fasta'}
```
- anno_outpath: output folder of the MegaFun annotation
- chocoPath: we use full_chocophlan.v201901_v31.tar.gz which can be downloaded from huttenhower.sph.harvard.edu, it is about 16GB
```
wget https://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
tar -xzf full_chocophlan.v201901_v31.tar.gz
```
- chocoMap, ncbi_name are prebuilt and provided in TaxChocoDB folder of MegaFun which are created by MegaFun/buildAnno/constructTaxChocoDB.R . 
- ncbiSummary in TaxChocoDB was downloaded from ncbi (September 2025) https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
-ncore: number of CPUs
-metaphlan_outpath: path to the output of Metaphlan
-metaphlan_prop=0.5: keep only species with abundance greater than metaphlan_prop


-runMetaPhlan: if you want to run MetaPhlan (YES) or not (NOT). If YES, then next parameters (sample_Dir, suffix_sample, file_type) are required. If you are familar with MetaPhlan, it is recommended to run MetaPhlan yourself in adavance.
-sample_Dir: the path to the input metagenomics data
-suffix_sample: need to specify the suffix: if it is compressed format (gz) or not
-file_type: parameter of Metaphlan: need to specify input format (fasta or fastq)


#### Of Note: generation of annotation include sequential steps which are wrapped up in MegaFun_buildAnno.sh. This step can take time and sometimes get issues of memory. In our experience, we was successful for all analyese with an allocation of 200GB memory. Experienced users should have a close look at those steps inside the following files: MegaFun_anno_pipeline.sh and MegaFun_buildAnno.sh, and might adapt to run them efficiently


### To run MegaFun for quantification:

```
bash /path/to/MegaFun/MegaFun.sh -in /path/to/yourinputdata -iso /path/to/isolate_genome_MegaFun/ -pan /path/to/pangenome_shortNames/ -smap /path/to/species_list_summary.txt -p 16
```

-species_list_summary.txt is the file containing the mapping between chocophlan and nbci which is created by buildAnno/extract_isolate_list.R. The file is already generated during the generation of MegaFun annotation. It is usually located at the output folder (anno_outpath: annoOut/species_list_summary.txt).


## output: MegaFun quantification produces an output folder:
- final_paralog.RData, final_singleton.RData for isolate level quantification
- Pan/XAEM_count for pangenome level quantification

## License  

GNU General Public License v3.0