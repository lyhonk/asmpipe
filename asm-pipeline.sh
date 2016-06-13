#!/bin/bash

usage()
{
	echo "usage: `basename $0` asm.config input.bam"
}

save_progress()
{
	progress=$1
	
	if [ -z "`grep -i error ${log_file}`" ] && [ -z "`grep -i "No such file or directory" ${log_file}`" ]
	then
		echo "${progress}=OK" >> ${progress_file} 		
	else
		echo "ERROR : current progress \"${progress}\""
		exit -1
	fi
}

check_tagmeth_asm_options()
{
	option=$2
	if [ -n "${!option}" ]
	then
		tagmeth_asm_options=${tagmeth_asm_options}" $1 ${!option}"
	fi
}

generate_trackhub()
{
	trackhub_genome_path=${trackhub_path}"/"${genome_type}
	mkdir -p ${trackhub_genome_path}
	trackhub_data_path=${trackhub_genome_path}"/data"
	mkdir -p ${trackhub_data_path}
	mv ${tagmeth_index_meth_bw} ${trackhub_data_path}"/" 
	mv ${tagmeth_index_umeth_bw} ${trackhub_data_path}"/"
	mv ${tagmeth_asm_bb} ${trackhub_data_path}"/"

	# hub.txt 

	hub_file=${trackhub_path}"/hub.txt"
	echo "hub hub_name(change this to your projectname)" > ${hub_file}
	echo "shortLabel hub_short_label" >> ${hub_file}
	echo "longLabel hub_long_label" >> ${hub_file}
	echo "genomesFile genomes.txt" >> ${hub_file}
	echo "email ${USER}_lilab@icsc.dlmedu.edu.cn" >> ${hub_file}
	echo "descriptionUrl hub.html" >> ${hub_file}

	# genomes.txt

	genomes_file=${trackhub_path}"/genomes.txt"
	echo "genome ${genome_type}" > ${genomes_file}
	echo "trackDb ${genome_type}/trackDb.txt" >> ${genomes_file}

	# hub.html

	hub_html=${trackhub_path}"/hub.html"
	echo "<!DOCTYPE html>" > ${hub_html}
	echo "<html>" >> ${hub_html}
	echo "<head>" >> ${hub_html}
	echo "Track Hub Description" >> ${hub_html}
	echo "</head>" >> ${hub_html}
	echo "<body>" >> ${hub_html}
	echo "<p> Your track hub description here. </p>" >> ${hub_html}
	echo "</body>" >> ${hub_html}
	echo "</html>" >> ${hub_html}
	
	# trackDb.txt
	trackdb_txt=${trackhub_genome_path}"/trackDb.txt"
	track_name="${input_bam_filebase}_TagMeth"
	echo "track ${track_name}" > ${trackdb_txt}
	echo "superTrack on show" >> ${trackdb_txt}
	echo "shortLabel ${track_name}" >> ${trackdb_txt}
	echo "longLabel ${track_name}" >> ${trackdb_txt}
	echo "html ${hub_html}" >> ${trackdb_txt}

	# subtracks
	echo "" >> ${trackdb_txt}
	echo "	track ${input_bam_filebase}_meth" >> ${trackdb_txt}
	echo "	parent ${track_name}" >> ${trackdb_txt}
	echo "	bigDataUrl	./data/`basename ${tagmeth_index_meth_bw}`" >> ${trackdb_txt}
	echo "	shortLabel	meth" >> ${trackdb_txt}
	echo "	longLabel	${input_bam_filebase}_meth" >> ${trackdb_txt}
	echo "	type bigWig" >> ${trackdb_txt}
	echo "	visibility hide" >> ${trackdb_txt}
	echo "	color 3,169,244" >> ${trackdb_txt}
	echo "	aggregate transparentOverlay" >> ${trackdb_txt}
	echo "	showSubtrackColorOnUi on" >> ${trackdb_txt}
	echo "	yLineOnOff on" >> ${trackdb_txt}
	echo "	priority 1.0" >> ${trackdb_txt}
	
	echo "" >> ${trackdb_txt}	
	echo "	track ${input_bam_filebase}_unmeth" >> ${trackdb_txt}
	echo "	parent ${track_name}" >> ${trackdb_txt}
	echo "	bigDataUrl	./data/`basename ${tagmeth_index_umeth_bw}`" >> ${trackdb_txt}
	echo "	shortLabel	unmeth" >> ${trackdb_txt}
	echo "	longLabel	${input_bam_filebase}_unmeth" >> ${trackdb_txt}
	echo "	type bigWig" >> ${trackdb_txt}
	echo "	visibility hide" >> ${trackdb_txt}
	echo "	color 139,195,74" >> ${trackdb_txt}
	echo "	aggregate transparentOverlay" >> ${trackdb_txt}
	echo "	showSubtrackColorOnUi on" >> ${trackdb_txt}
	echo "	yLineOnOff on" >> ${trackdb_txt}
	echo "	priority 1.0" >> ${trackdb_txt}

	echo "" >> ${trackdb_txt}
	echo "	track ${input_bam_filebase}_asm" >> ${trackdb_txt}
	echo "	parent ${track_name}" >> ${trackdb_txt}
	echo "	bigDataUrl	./data/`basename ${tagmeth_asm_bb}`" >> ${trackdb_txt}
	echo "	shortLabel	asm" >> ${trackdb_txt}
	echo "	longLabel	${input_bam_filebase}_asm" >> ${trackdb_txt}
	echo "	type bigBed" >> ${trackdb_txt}
	echo "	visibility hide" >> ${trackdb_txt}
	echo "	color 255,193,7" >> ${trackdb_txt}
	echo "	aggregate transparentOverlay" >> ${trackdb_txt}
	echo "	showSubtrackColorOnUi on" >> ${trackdb_txt}
	echo "	yLineOnOff on" >> ${trackdb_txt}
	echo "	priority 1.0" >> ${trackdb_txt}
}

if [ $# != 2  ] 
then
	usage
	exit -1
fi

asm_config_file=$1
if [ ! -f ${asm_config_file} ]
then
	echo "Error: \"${asm_config_file}\" must be a valid asm config file."
	exit -1
fi
source ${asm_config_file}

input_bam_filename=$2
if [ ! -f ${input_bam_filename} ]
then
	echo "Error: \"${input_bam_filename}\" must be a valid bam file."
	exit -1
fi

# tool cmds

script_path=`dirname $0`
chromsize_path="${script_path}/chrom.sizes"
python_cmd="python"
rmdump_cmd="${python_cmd} ${script_path}/rmdup.py -s -q"
tagmeth_cmd="${python_cmd} ${script_path}/tagmeth.py -q"
asm_cmd="${script_path}/tagmeth.R"
wigtobigwig_cmd="${script_path}/utils/wigToBigWig"
bedtobigbed_cmd="${script_path}/utils/bedToBigBed"

# file path & names

reffa_file=$2
input_bam_filebase=`basename ${input_bam_filename}`
input_bam_filebase=${input_bam_filebase%.*}
log_file="./${input_bam_filebase}.pipline.log"
progress_file=".${input_bam_filebase}.progress"

# read config file

source ${asm_config_file}

if [ ! -f ${ref_fa} ]
then 
	echo " Invalid reference file ref_fa=${ref_fa} in ${conf_file} "
	exit -1
fi

if [ "${genome_type}" != "mm9"  ] && [ "${genome_type}" != "mm10" ] && [ "${genome_type}" != "hg18" ] && [ "${genome_type}" != "hg19" ] 
then
	echo "    Error: unknown genome_type - ${genome_type}, only mm9 | mm10 | hg18 | hg19 are surported"
	exit -1
fi

# Initialize
echo `date` | tee ${log_file}
echo "[*] Preparing output files & path..." | tee -a ${log_file}

# prepare progress file

if [ -f ${progress_file} ]
then
	source ${progress_file}
else
	touch ${progress_file}
fi

# output path

output_path="${input_bam_filebase}.tagmeth"
if [ ! -d ${output_path} ]
then
	mkdir ${output_path}
fi

trackhub_path="${input_bam_filebase}.trackhub"
if [ ! -d ${trackhub_path} ]
then
	mkdir ${trackhub_path}
fi

# step 1: clean up reads

#echo "[*] Cleaning up reads..." | tee -a ${log_file}
#current_progress="asm_clean_reads"
#clean_bam_filename="./${output_path}/clean.${input_bam_filebase}.bam"

#if [ -n  ${key_re} ]
#then
#	rmdump_cmd=${rmdump_cmd}" -k ${key_re} "
#fi

#if [ -n ${num_re} ]
#then
#	rmdump_cmd=${rmdump_cmd}" -n ${num_re} "
#fi

#if [ -z ${!current_progress} ] 
#then
#	cmd="${rmdump_cmd} -o ${clean_bam_filename} ${input_bam_filename}"
#	echo "    CMD - ${cmd}" | tee -a ${log_file}
#	${cmd} >> ${log_file} 2>&1
#
#	save_progress ${current_progress}
#fi

# step 2: generate tag meth csv file
echo `date`| tee -a ${log_file}
echo "[*] Generating tagmeth..." | tee -a ${log_file}
current_progress="asm_generate_tagmeth"
tagmeth_file="./${output_path}/${input_bam_filebase}.tagmeth.csv"
if [ -z ${!current_progress} ]
then
	cmd="${tagmeth_cmd} ${input_bam_filename} ${ref_fa} -o ${tagmeth_file}"
	echo "    CMD - ${cmd}" | tee -a ${log_file}
	${cmd} >> ${log_file} 2>&1
	
	save_progress ${current_progress}	
fi

# step 3: cleanup & split tagmeth by chrom

echo "[*] Cleaning up & spliting tagmeth" | tee -a ${log_file}
current_progress="asm_clean_tagmeth"
#clean_tagmeth_file="./${output_path}/clean.${input_bam_filebase}.tagmeth.csv"
chr_tagmeth_path="./${output_path}/${input_bam_filebase}.csv"
if [ -z ${!current_progress} ]
then
	# clean up tagmeth file

#	keysize=`tail -1 ${tagmeth_file} | cut -f1 -d" " | wc -c`
#	echo "    removing duplicated tagmeth" | tee -a ${log_file}
#	uniq -w ${keysize} ${tagmeth_file} > ${clean_tagmeth_file}
	
	# split tagmeth by chrom

	if [ ! -d ${chr_tagmeth_path} ]
	then	
		mkdir -p ${chr_tagmeth_path}
	fi	
	rm -rf ${chr_tagmeth_path}/*
	
	awk -v path=${chr_tagmeth_path} ' NR!=1 {print >> path"/"$4".tagmeth.csv"}' ${tagmeth_file} 

	save_progress ${current_progress}
fi

# step 4: 
# a) generate tagmeth index 
# b) find asm using ALE method 
# c) generate figures & trackhubs files
echo `date` | tee -a ${log_file}
echo "[*] Processing tagmeth files" | tee -a ${log_file}
current_progress="asm_process_tagmeth"

tagmeth_asm_options="-t ${genome_type}"
check_tagmeth_asm_options "-m" "meth_score"
check_tagmeth_asm_options "-u" "unmeth_score"
check_tagmeth_asm_options "-w" "sw_size"
check_tagmeth_asm_options "-a" "asm_size"
check_tagmeth_asm_options "-i" "meth_index"
check_tagmeth_asm_options "-x" "unmeth_index"

if [ -z ${!current_progress} ]
then
	cmd="${asm_cmd} ${tagmeth_asm_options} ${chr_tagmeth_path}"	
	echo "    CMD - ${cmd}" | tee -a ${log_file}
	${cmd} >> ${log_file} 2>&1
	
	# save rds files		

	chr_tagmeth_rds_path="./${output_path}/${input_bam_filebase}.rds"
	rm -rf ${chr_tagmeth_rds_path}
	mkdir -p ${chr_tagmeth_rds_path}	
	mv ./*.rds "${chr_tagmeth_rds_path}/"	
	mv ./*.csv "./${output_path}/"	
	
	save_progress ${current_progress}
fi

# generate trackhub files
echo `date` | tee -a ${log_file}
echo "[*] Generating trackhub files" | tee -a ${log_file}
current_progress="generate_trackhub"
if [ -z ${!current_progress} ]
then
	chrom_sizes="${chromsize_path}/${genome_type}.chrom.sizes"
	tagmeth_index_meth_wig="./${input_bam_filebase}.idx.meth.wig"
	tagmeth_index_umeth_wig="./${input_bam_filebase}.idx.umeth.wig"
	tagmeth_asm_bed="./${input_bam_filebase}.asm.bed"
	
	tagmeth_index_meth_bw="./${input_bam_filebase}.idx.meth.bw"
	tagmeth_index_umeth_bw="./${input_bam_filebase}.idx.umeth.bw"
	tagmeth_asm_bb="./${input_bam_filebase}.asm.bb"	

	${wigtobigwig_cmd} ${tagmeth_index_meth_wig} ${chrom_sizes} ${tagmeth_index_meth_bw}  >> ${log_file} 2>&1  
	${wigtobigwig_cmd} ${tagmeth_index_umeth_wig} ${chrom_sizes} ${tagmeth_index_umeth_bw}  >> ${log_file} 2>&1  
	
	tagmeth_sorted_asm_bed="./sorted.${input_bam_filebase}.asm.bed"
	sort -k1,1 -k2,2n ${tagmeth_asm_bed} > ${tagmeth_sorted_asm_bed}
	${bedtobigbed_cmd} ${tagmeth_sorted_asm_bed} ${chrom_sizes} ${tagmeth_asm_bb}
	rm ${tagmeth_sorted_asm_bed}	

	generate_trackhub

	save_progress ${current_progress}
fi

# final

mv *.wig *.bed ${output_path}
echo `date` | tee -a ${log_file}
echo "[*] Complete `date`" | tee -a ${log_file} 


