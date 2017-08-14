#!/usr/bin/env nextflow
// CIPHER Main Script
// Version 1.0.0
// Author: Carlos Guzman
// E-mail: cag104@ucsd.edu

// default modules
params.help = false
params.subsample = false

// default parameters
params.threads = 1
params.outdir = './results'
params.subsampled_reads = 100000

// default clumpify parameters
params.clumpify = true

// default bbduk parameters
params.bbduk = true
params.bbduk_ktrim = 'r'
params.bbduk_k = 23
params.bbduk_mink = 11
params.bbduk_hdist = 1
params.bbduk_trimq = 6
params.bbduk_minlength = 10

// default fastqc parameters
params.fastqc = true

// default multiqc parameters
params.multiqc = true

// default dge parameters
params.dge = false

// default mapping parameters
params.mapping = true
params.aligner = 'bbmap'

params.bbmap_maxindel = 1
params.bbmap_ambig = "random"
params.bbmap_minid = 0.76
params.bbmap_intronlen = 20

params.bwa_T = 30
params.bwa_k = 19
params.bwa_w = 100
params.bwa_d = 100
params.bwa_r = 1.5
params.bwa_c = 10000
params.bwa_A = 1
params.bwa_B = 4
params.bwa_O = 6
params.bwa_E = 1
params.bwa_L = 5
params.bwa_U = 9

params.bt2_D = 20
params.bt2_R = 3
params.bt2_N = 0
params.bt2_L = 20
params.bt2_i = 'S,5,1,0.50'
params.bt2_trim5 = 0
params.bt2_trim3 = 0

params.hs_trim5 = 0
params.hs_trim3 = 0
params.hs_mp = '6,2'
params.hs_sp = '2,1'
params.hs_np = 1
params.hs_rdg = '5,3'
params.hs_rfg = '5,3'
params.hs_pen_cansplice = 0
params.hs_pen_noncansplice = 12
params.hs_min_intronlen = 20
params.hs_max_intronlen = 500000
params.hs_k = 5
params.hs_max_seeds = 5

params.star_clip3pNbases = 0
params.star_clip5pNbases = 0
params.star_outFilterMultimapScoreRange = 1
params.star_outFilterMultimapNmax = 10
params.star_outFilterMismatchNmax = 10
params.star_outFilterScoreMin = 0
params.star_alignEndsType = 'Local'
params.star_winAnchorMultimapNmax = 50
params.star_quantMode = '-'
params.star_twopassMode = 'None'

// default bamcoverage parameters
params.bamcoverage = true
params.bamcoverage_bs = 10
params.bamcoverage_smooth = 100
params.bamcoverage_e = 150

// default epic parameters
params.epic_w = 200
params.epic_g = 3
params.epic_qvalue = 0.01

// default macs parameters
params.macs_qvalue = 0.01

// default danpos2 parameters
params.danpos_jw = 40
params.danpos_jd = 150

// default downstream analysis parameters
params.downstream_analysis = true



// print help
if (params.help == true) {
	log.info ''
	log.info '========================================================================================'
	log.info ''
	log.info 'Usage:'
	log.info 'nextflow run cipher_main.nf --mode <MODE> --config <CONFIG> --fa <FASTA> --gtf <GTF> --lib <LIB> --readLen <LEN> [options]'
	log.info ''
	log.info 'REQUIRED FLAGS:'
	log.info ''
	log.info '--mode				Choose from available: chip, rna, gro, mnase, dnase, atac, analysis.'
	log.info '--config			A tab separated text file containing sample information. Check README for more information.'
	log.info '--fa				A reference genome in FASTA format.'
	log.info '--gtf				Reference genome in GTF format.'
 	log.info '--lib				Library information. "s" for single-stranded data, "p" for pair-ended data.'
 	log.info '--readLen			The length of your reads.'
	log.info ''
	log.info 'CLUMPIFY FLAGS:'
	log.info ''
	log.info '--clumpify			Set this to false if you would like to skip the Clumpify compression step. (Default: true)'
	log.info ''
	log.info 'BBDUK FLAGS:'
	log.info ''
	log.info '--bbduk				Set this to false if you would like to skip the BBDuk adapter trimming step. (Default: true)'
	log.info ''
	log.info '--bbduk_ktrim			See ktrim flag in BBMap user manual for more information. (Default: r)'
	log.info '--bbduk_k				See k flag in BBMap user manual for more information. (Default: 23)'
	log.info '--bbduk_mink			See mink flag in BBMap user manual for more information. (Default: 11)'
	log.info '--bbduk_hdist			See hdist flag in BBMap user manual for more information. (Default: 1)'
	log.info '--bbduk_trimq			See trimq flag in BBMap user manual for more information. (Default: trimq)'
	log.info '--bbduk_minlength		See minlength flag in BBMap user manual for more information. (Default: 10)'
	log.info ''
	log.info 'FASTQC FLAGS:'
	log.info ''
	log.info '--fastqc			Set this to false if you would like to skip the FastQC step. (Default: true)'
	log.info ''
	log.info 'MULTIQC FLAGS:'
	log.info ''
	log.info '--multiqc			Set this to false if you would like to skip the MultiQC step. (Default: true)'
	log.info ''
	log.info 'DGE FLAGS:'
	log.info ''
	log.info '--dge				Set this to false if you would like to skip the differential gene expression step. (Default: false)'
	log.info ''
	log.info 'ALIGNMENT FLAGS:'
	log.info ''
	log.info '--mapping			Set this to false if you would like to skip the alignment step. (Default: true)'
	log.info '--aligner			The aligner to map your data to the reference genome. Choose from: bbmap, bowtie2, bwa, hisat2, star. (Default: bbmap)'
	log.info ''
	log.info '--bbmap_maxindel	See maxindel flag in BBMap user manual for more information. Set to 200k for RNA-seq. (Default: 1)'
	log.info ''
	log.info '--bwa_T			See T flag in BWA user manual for more information. (Default: 30)'
	log.info '--bwa_k			See k flag in BWA user manual for more information. (Default: 19)'
	log.info '--bwa_w			See w flag in BWA user manual for more information. (Default: 100)'
	log.info '--bwa_d			See d flag in BWA user manual for more information. (Default: 100'
	log.info '--bwa_r			See r flag in BWA user manual for more information. (Default: 1.5)'
	log.info '--bwa_c			See c flag in BWA user manual for more information. (Default: 10000)'
	log.info '--bwa_A			See A flag in BWA user manual for more information. (Default: 1)'
	log.info '--bwa_B			See B flag in BWA user manual for more information. (Default: 4)'
	log.info '--bwa_O			See O flag in BWA user manual for more information. (Default: 6)'
	log.info '--bwa_E			See E flag in BWA user manual for more information. (Default: 1)'
	log.info '--bwa_L			See L flag in BWA user manual for more information. (Default: 5)'
	log.info '--bwa_U			See U flag in BWA user manual for more information. (Default: 9)'
	log.info ''
	log.info '--bt2_D			See D flag in Bowtie2 user manual for more information. (Default: 20)'
	log.info '--bt2_R			See R flag in Bowtie2 user manual for more information. (Default: 3)'
	log.info '--bt2_N			See N flag in Bowtie2 user manual for more information. (Default: 0)'
	log.info '--bt2_L			See L flag in Bowtie2 user manual for more information. (Default: 20)'
	log.info '--bt2_i			See i flag in Bowtie2 user manual for more information. (Default: S,5,1,0.50)'
	log.info '--bt2_trim5			See trim5 flag in Bowtie2 user manual for more information. (Default: 0)'
	log.info '--bt2_trim3			See trim3 flag in Bowtie2 user manual for more information. (Default: 0)'
	log.info ''
	log.info 'BAMCOVERAGE FLAGS:'
	log.info ''
	log.info '--bamcoverage			Set this to false if you would like to skip the bamCoverage bigwig creation step. (Default: true)'
	log.info ''
	log.info 'OTHER FLAGS:'
	log.info ''
	log.info '--threads			The number of threads run per process. (Default: 1)'
	log.info '--outdir			The name of the output directory. (Default: ./results)'
	log.info ''
	log.info 'CHIP-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info '--macs_g					See g flag in MACS2 user manual for more information. (Default: automatically calculated)'
	log.info '--epic_egs				See egs flag in EPIC user manual for more information. (Default: automatically calculated)'
	log.info ''
	log.info 'DNASE-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info 'ATAC-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info 'MNASE-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info 'RNA-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info 'GRO-SEQ FLAGS:'
	log.info ''
	log.info '--downstream_analysis		Set this to false if you would like CIPHER to skip the downstream analysis. (Default: true)'
	log.info ''
	log.info 'TESTING FLAGS:'
	log.info ''
	log.info '--subsample			Set this to true if you would like to subsample raw fastq files into smaller subsets for quick pipeline testing. (Default: false)'
	log.info '--subsample_reads		The number of sampled reads. (Default: 100k)'
	log.info ""
	log.info '========================================================================================'
	exit 1
}

// parameter check
 if (!params.config) {
 	exit 1, "Please specify a config file."
 }

 if (!params.fa) {
 	exit 1, "Please specify a FASTA file."
 }

 if (!params.gtf) {
 	exit 1, "Please specify a GTF file."
 }

 if (!params.mode) {
 	exit 1, "Please specify a mode. Choose from: 'chip', 'rna', 'gro', 'dnase', 'mnase', 'atac', or 'analysis'."
 }

 if (!params.lib) {
 	exit 1, "Please specify the type of sequencing library your config file contains. Available: 's' for single-end sequencing and 'p' for pair-end sequencing."
 }

 if (!params.readLen) {
 	exit 1, "Please specify the length of your reads."
 }

 // set up file shortcuts
 config_file = file(params.config)
 fasta_file = file(params.fa)
 gtf_file = file(params.gtf)

 // display run information
log.info ""
log.info "                 RUN INFORMATION               "
log.info "==============================================="
log.info ''
log.info "mode:				${params.mode}"
log.info "config:				${config_file}"
log.info "fasta:				${fasta_file}"
log.info "gtf:				${gtf_file}"
log.info "library:			${params.lib}"
log.info "read-length:			${params.readLen}"
log.info ""
log.info "subsampled:			${params.subsample}"
log.info ""
log.info "clumpify:			${params.clumpify}"
log.info "bbduk:				${params.bbduk}"
log.info "fastqc:				${params.fastqc}"
log.info "multiqc:			${params.multiqc}"
log.info "mapping:			${params.mapping}"
log.info "bamcoverage:			${params.bamcoverage}"
log.info "downstream_analysis		${params.downstream_analysis}"
log.info ""
log.info "aligner:			${params.aligner}"
log.info "==============================================="
log.info ""

// Parse single ended config file
if (params.lib == "s" || params.lib == "S") {
fastqs_s = Channel
 	.from(config_file.readLines())
 	.map { line ->
 			list = line.split()
 			mergeid = list[0]
 			id = list[1]
 			path1 = file(list[2])
 			controlid = list[3]
 			mark = list[4]
 			[ mergeid, id, path1, controlid, mark ]
 		}
 	}

// Parse pair ended config file
if (params.lib == "p" || params.lib == "P") {
fastqs_p = Channel
 	.from(config_file.readLines())
 	.map { line ->
 			list = line.split()
 			mergeid = list[0]
 			id = list[1]
 			path1 = file(list[2])
 			path2 = file(list[3])
 			controlid = list[4]
 			mark = list[5]
 			[ mergeid, id, path1, path2, controlid, mark ]
 		}
 	}

// subsample data if required for single end data
if (params.subsample == true && params.lib == "s") {
	process subsample {

		input:
		set mergeid, id, file(read1), controlid, mark from fastqs_s

		output:
		set mergeid, id, file("${id}.subsampled.fq.gz"), controlid, mark into subsampled_fastqs_s

		script:
		"""
		reformat.sh in=${read1} out=${id}.subsampled.fq.gz sample=${params.subsampled_reads}
		"""
	}

	subsampled_fastqs_s.set {

		clumpify_fqs_s
	}
} else if (params.subsample == false && params.lib == "s") {
	fastqs_s.set {

		clumpify_fqs_s
	}
}

// subsample data if required for pair end data
if (params.subsample == true && params.lib == "p") {
	process subsample {

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_p

		output:
		set mergeid, id, file("${id}_R1.subsampled.fq.gz"), file("${id}_R2.subsampled.fq.gz"), controlid, mark into subsampled_fastqs_p

		script:
		"""
		reformat.sh in=${read1} in2=${read2} out=${id}_R1.subsampled.fq.gz out2=${id}_R2.subsampled.fq.gz sample=${params.subsampled_reads}
		"""
	}

	subsampled_fastqs_p.set {

		clumpify_fqs_p
	}
} else if (params.subsample == false && params.lib == "p") {
	fastqs_p.set {

		clumpify_fqs_p
	}
}

// Fetch chromosome sizes from FASTA
if (params.downstream_analysis == true)
    process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic, chrom_sizes_danpos

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

// Calculate effective genome size
if (!params.macs_g && !params.epic_egs && params.downstream_analysis == true) {
	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file
 		file "egs_ratio.txt" into egs_ratio

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_macs_input
 			egs_size_macs_no_input
 		}

egs_ratio.map{ file ->
 			file.text.trim() } .set {
 				egs_ratio
 			}

 			egs_ratio.set {
 				egs_ratio_epic
 			}
 		}

// create indexes for different aligners
if (params.aligner == 'bbmap' && params.mapping == true) {
 	   // Generate BBMap Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file}
 		"""
 	}
 }

if (params.aligner == 'bowtie2' && params.mapping == true) {
	// Generate Bowtie2 Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/bowtie2_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("*") into bowtie2_index

 		script:
 		"""
 		bowtie2-build --threads ${params.threads} ${fasta_file} genome
 		"""
 	}
 }

if (params.aligner == 'bwa' && params.mapping == true) {
	// Generate BWA Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/bwa_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("*") into bwa_index

 		script:
 		"""
 		bwa index -p genome ${fasta_file}
 		"""
 	}
 }

if (params.aligner == 'star' && params.mapping == true) {
	// Generate STAR Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/star_index", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file

 		output:
 		file("indexFiles/*") into star_index

 		script:
 		overhang = Math.round((${params.readLen} as int) - 1)
 		"""
 		mkdir indexFiles
 		STAR --runMode genomeGenerate --runThreadN ${params.threads} --genomeDir indexFiles --genomeFastaFiles ${fasta_file} --sjdbGTFfile ${gtf_file} --sjdbOverhang ${overhang}
 		"""
 	}
 }

if (params.aligner == 'hisat2' && params.mapping == true) {
	// Generate HISAT2 Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/hisat2_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("*") into hisat2_index

 		script:
 		"""
 		hisat2-build -p ${params.threads} ${fasta_file} genome
 		"""
 	}
 }

// step clumpify for single end data
if (params.clumpify == true && params.lib == "s") {
	process clumpify {

	input:
	set mergeid, id, file(read1), controlid, mark from clumpify_fqs_s

	output:
	set mergeid, id, file("${id}.clumped.fq.gz"), controlid, mark into pre_fastqc_fqs_s, bbduk_fqs_s
	file("clumpify_parameters_${id}.txt")

	script:
	"""
	clumpify.sh -Xmx12g in=${read1} out=${id}.clumped.fq.gz
	echo 'clumpify.sh -Xmx12g in=${read1} out=${id}.clumped.fq.gz' > clumpify_parameters_${id}.txt
	"""
}

	} else if (params.clumpify == false && params.lib == "s") {
		clumpify_fqs_s.into {

			pre_fastqc_fqs_s
			bbduk_fqs_s
		}
	}

// step clumpify for pair end data
if (params.clumpify == true && params.lib == "p") {
	process clumpify {

	input:
	set mergeid, id, file(read1), file(read2), controlid, mark from clumpify_fqs_p

	output:
	set mergeid, id, file("${id}_R1.clumped.fq.gz"), file("${id}_R2.clumped.fq.gz"), controlid, mark into pre_fastqc_fqs_p, bbduk_fqs_p
	file("clumpify_parameters_${id}.txt")

	script:
	"""
	clumpify.sh -Xmx12g in=${read1} in2=${read2} out=${id}_R1.clumped.fq.gz out2=${id}_R2.clumped.fq.gz
	echo 'clumpify.sh -Xmx12g in=${read1} in2=${read2} out=${id}_R1.clumped.fq.gz out2=${id}_R2.clumped.fq.gz' > clumpify_parameters_${id}.txt
	"""
}

	} else if (params.clumpify == false && params.lib == "p") {
		clumpify_fqs_p.into {

			pre_fastqc_fqs_p
			bbduk_fqs_p
		}
	}

// step bbduk for single end data
if (params.bbduk == true && params.lib == "s") {
	process bbduk {

		publishDir "${params.outdir}/${params.mode}/${id}/trimmed_reads", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from bbduk_fqs_s

		output:
		set mergeid, id, file("${id}.trimmed.fq.gz"), controlid, mark into post_fastqc_fqs_s, mapping_fqs
		file("${id}_bbduk_statsFile.txt")
		file("bbduk_parameters_${id}.txt")

		script:
		"""
		bbduk.sh in=${read1} out=${id}.trimmed.fq.gz ref=$baseDir/adapters/adapters.fa ktrim=${params.bbduk_ktrim} k=${params.bbduk_k} mink=${params.bbduk_mink} hdist=${params.bbduk_hdist} trimq=${params.bbduk_trimq} minlength=${params.bbduk_minlength} tpe ordered=t stats=${id}_bbduk_statsFile.txt
		echo 'bbduk.sh in=${read1} out=${id}.trimmed.fq.gz ref=$baseDir/adapters/adapters.fa ktrim=${params.bbduk_ktrim} k=${params.bbduk_k} mink=${params.bbduk_mink} hdist=${params.bbduk_hdist} trimq=${params.bbduk_trimq} minlength=${params.bbduk_minlength} tpe ordered=t stats=${id}_bbduk_statsFile.txt' > bbduk_parameters_${id}.txt
		"""
	}
} else if (params.bbduk == false && params.lib == "s") {
	bbduk_fqs_s.into {

		post_fastqc_fqs_s
		mapping_fqs
	}
}

// step bbduk for pair end data
if (params.bbduk == true && params.lib == "p") {
	process bbduk {

		publishDir "${params.outdir}/${params.mode}/${id}/trimmed_reads", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from bbduk_fqs_p

		output:
		set mergeid, id, file("${id}_R1.trimmed.fq.gz"), file("${id}_R2.trimmed.fq.gz"), controlid, mark into post_fastqc_fqs_p, mapping_fqs
		file("${id}_bbduk_statsFile.txt")
		file("bbduk_parameters_${id}.txt")

		script:
		"""
		bbduk.sh in=${read1} in2=${read2} out=${id}_R1.trimmed.fq.gz out2=${id}_R2.trimmed.fq.gz ref=$baseDir/adapters/adapters.fa ktrim=${params.bbduk_ktrim} k=${params.bbduk_k} mink=${params.bbduk_mink} hdist=${params.bbduk_hdist} trimq=${params.bbduk_trimq} minlength=${params.bbduk_minlength} tpe tbo ordered=t stats=${id}_bbduk_statsFile.txt
		echo 'bbduk.sh in=${read1} in2=${read2} out=${id}_R1.trimmed.fq.gz out2=${id}_R2.trimmed.fq.gz ref=$baseDir/adapters/adapters.fa ktrim=${params.bbduk_ktrim} k=${params.bbduk_k} mink=${params.bbduk_mink} hdist=${params.bbduk_hdist} trimq=${params.bbduk_trimq} minlength=${params.bbduk_minlength} tpe tbo ordered=t stats=${id}_bbduk_statsFile.txt' > bbduk_parameters_${id}.txt
		"""
	}
} else if (params.bbduk == false && params.lib == "p") {
	bbduk_fqs_p.into {

		post_fastqc_fqs_p
		mapping_fqs
	}
}

// step pre fastqc for single end data
if (params.fastqc == true && params.lib == "s") {
	process pre_fastqc {

		publishDir "${params.outdir}/${params.mode}/${id}/fastqc", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from pre_fastqc_fqs_s

		output:
		file("*.{zip, html}") into pre_fastqc_multiqc_s

		script:
		"""
		fastqc ${read1}
		"""
	}
}

// step pre fastqc for pair end data
if (params.fastqc == true && params.lib == "p") {
	process pre_fastqc {

		publishDir "${params.outdir}/${params.mode}/${id}/fastqc", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from pre_fastqc_fqs_p

		output:
		file("*.{zip, html}") into pre_fastqc_multiqc_p

		script:
		"""
		fastqc -t 2 ${read1} ${read2}
		"""
	}
}

// step post fastqc for single end data
if (params.fastqc == true && params.lib == "s" && params.bbduk == true) {
	process post_fastqc {

		publishDir "${params.outdir}/${params.mode}/${id}/fastqc", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from post_fastqc_fqs_s

		output:
		file("*.{zip, html}") into post_fastqc_multiqc_s

		script:
		"""
		fastqc ${read1}
		"""
	}
}

// step post fastqc for pair end data
if (params.fastqc == true && params.lib == "p" && params.bbduk == true) {
	process post_fastqc {

		publishDir "${params.outdir}/${params.mode}/${id}/fastqc", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from post_fastqc_fqs_p

		output:
		file("*.{zip, html}") into post_fastqc_multiqc_p

		script:
		"""
		fastqc -t 2 ${read1} ${read2}
		"""
	}
}

// step mapping using bbmap single end data
if (params.mapping == true && params.aligner == "bbmap" && params.lib == "s") {
	process mapping_bbmap {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from mapping_fqs
		file("ref/*") from bbmap_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("${id}.bbmap_alignmentReport.txt")
		file("${id}.unmapped.bam")
		file("bbmap_parameters_${id}.txt")

		script:
		"""
		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam maxindel=${params.bbmap_maxindel} intronlen=${params.bbmap_intronlen} ambig=${params.bbmap_ambig} statsfile=${id}.bbmap_alignmentReport.txt minid=${params.bbmap_minid}
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam maxindel=${params.bbmap_maxindel} intronlen=${params.bbmap_intronlen} ambig=${params.bbmap_ambig} statsfile=${id}.bbmap_alignmentReport.txt minid=${params.bbmap_minid}' > bbmap_parameters_${id}.txt
		"""
	}
}

// step mapping using bbmap pair end data
if (params.mapping == true && params.aligner == "bbmap" && params.lib == "p") {
	process mapping_bbmap {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from mapping_fqs
		file("ref/*") from bbmap_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("${id}.bbmap_alignmentReport.txt")
		file("${id}.unmapped.bam")
		file("bbmap_parameters_${id}.txt")

		script:
		"""
		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam maxindel=${params.bbmap_maxindel} intronlen=${params.bbmap_intronlen} ambig=${params.bbmap_ambig} statsfile=${id}.bbmap_alignmentReport.txt minid=${params.bbmap_minid}
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bbmap.sh 'in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam maxindel=${params.bbmap_maxindel} intronlen=${params.bbmap_intronlen} ambig=${params.bbmap_ambig} statsfile=${id}.bbmap_alignmentReport.txt minid=${params.bbmap_minid}' > bbmap_parameters_${id}.txt
		"""
	}
}

// step mapping using bowtie2 single end data
if (params.mapping == true && params.aligner == "bowtie2" && params.lib == "s") {
	process mapping_bowtie2 {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from mapping_fqs
		file("*") from bowtie2_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("bowtie2_parameters_${id}.txt")
		file("${id}.bowtie2_alignmentReport.txt")
		file("bowtie2_parameters_${id}.txt")

		script:
		"""
		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} --trim5 ${params.bt2_trim5} --trim3 ${params.bt2_trim3} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam 2> ${id}.bowtie2_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} --trim5 ${params.bt2_trim5} --trim3 ${params.bt2_trim3} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam' > bowtie2_parameters_${id}.txt
		"""
	}
}

// step mapping using bowtie2 pair end data
if (params.mapping == true && params.aligner == "bowtie2" && params.lib == "p") {
	process mapping_bowtie2 {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from mapping_fqs
		file("*") from bowtie2_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("bowtie2_parameters_${id}.txt")
		file("${id}.bowtie2_alignmentReport.txt")
		file("bowtie2_parameters_${id}.txt")

		script:
		"""
		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} --trim5 ${params.bt2_trim5} --trim3 ${params.bt2_trim3} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam 2> ${id}.bowtie2_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} --trim5 ${params.bt2_trim5} --trim3 ${params.bt2_trim3} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam' > bowtie2_parameters_${id}.txt
		"""
	}
}

// step mapping using bwa single end data
if (params.mapping == true && params.aligner == "bwa" && params.lib == "s") {
	process mapping_bwa {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from mapping_fqs
		file("*") from bwa_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("bwa_parameters_${id}.txt")
		file("${id}.bwa_alignmentReport.txt")
		file("bwa_parameters_${id}.txt")

		script:
		"""
		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam 2> ${id}.bwa_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1}' > bwa_parameters_${id}.txt
		"""
	}
}

// step mapping using bwa pair end data
if (params.mapping == true && params.aligner == "bwa" && params.lib == "p") {
	process mapping_bwa {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from mapping_fqs
		file("*") from bwa_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, spp_bams, bam_grouping, bam_danpos, bam_qorts, bam_featurecounts, bam_stringtie
		file("bwa_parameters_${id}.txt")
		file("${id}.bwa_alignmentReport.txt")
		file("bwa_parameters_${id}.txt")

		script:
		"""
		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam 2> ${id}.bwa_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2}' > bwa_parameters_${id}.txt
		"""
	}
}

// step mapping using hisat2 single end data
if (params.mapping == true && params.aligner == "hisat2" && params.lib == "s") {
	process mapping_hisat2 {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from mapping_fqs
		file("*") from hisat2_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, bam_grouping, bam_danpos, bam_preseq, bam_qorts, bam_featurecounts, bam_stringtie
		file("hisat2_parameters_${id}.txt")
		file("${id}.hisat2_alignmentReport.txt")
		file("hisat2_parameters_${id}.txt")

		script:
		"""
		hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam 2> ${id}.hisat2_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam' > hisat2_parameters_${id}.txt
		"""
	}
}

// step mapping using hisat2 pair end data
if (params.mapping == true && params.aligner == "hisat2" && params.lib == "p") {
	process mapping_hisat2 {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from mapping_fqs
		file("*") from hisat2_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, bam_grouping, bam_danpos, bam_preseq, bam_featurecounts, bam_qorts, bam_stringtie
		file("hisat2_parameters_${id}.txt")
		file("${id}.hisat2_alignmentReport.txt")
		file("hisat2_parameters_${id}.txt")

		script:
		"""
		hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam 2> ${id}.hisat2_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam' > hisat2_parameters_${id}.txt
		"""
	}
}

// step mapping using star single end data
if (params.mapping == true && params.aligner == "star" && params.lib == "s") {
	process mapping_star {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), controlid, mark from mapping_fqs
		file("indexFiles/*") from star_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, bam_grouping, bam_danpos, bam_preseq, bam_featurecounts, bam_qorts, bam_stringtie
		file("star_parameters_${id}.txt")
		file("${id}.star_alignmentReport.txt")
		file("star_parameters_${id}.txt")

		script:
		"""
		STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode} 2> ${id}.star_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}Aligned.sortedByCoord.out.bam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode}' > star_parameters_${id}.txt
		"""
	}
}

// step mapping using star pair end data
if (params.mapping == true && params.aligner == "star" && params.lib == "p") {
	process mapping_star {

		publishDir "${params.outdir}/${params.mode}/${id}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(read1), file(read2), controlid, mark from mapping_fqs
		file("indexFiles/*") from star_index

		output:
		set mergeid, id, file("*.sorted.mapped.bam"), controlid, mark, file("*.bam.bai") into bamcoverage_bams, bam_grouping, bam_danpos, bam_preseq, bam_featurecounts, bam_qorts, bam_stringtie
		file("star_parameters_${id}.txt")
		file("${id}.star_alignmentReport.txt")
		file("star_parameters_${id}.txt")

		script:
		"""
		STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} ${read2} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode} 2> ${id}.star_alignmentReport.txt
		samtools sort -T ${id} -@ ${params.threads} -o ${id}.sorted.mapped.bam ${id}Aligned.sortedByCoord.out.bam
		samtools index -@ ${params.threads} ${id}.sorted.mapped.bam
		echo 'STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} ${read2} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode}' > star_parameters_${id}.txt
		"""
	}
}

// merge replicates together then create bigwigs
if (params.bamcoverage == true) {

	singleBams = Channel.create()
	groupedBams = Channel.create()

	bamcoverage_bams.groupTuple(by: [0,3,4])
	.choice(singleBams, groupedBams) {
		it[2].size() > 1 ? 1 : 0
	}

	process merge_bams {

		publishDir "${params.outdir}/${params.mode}/merged/${mergeid}/alignments", mode: 'copy'

		input:
		set mergeid, id, file(bam), controlid, mark, file(bam_index) from groupedBams

		output:
		set mergeid, id, file("${mergeid}_merged.bam"), controlid, mark, file("${mergeid}_merged.bam.bai") into mergedBams
		file("samtools_merge_parameters_${mergeid}.txt")

		script:
		def id = id.sort().join(':')
		"""
		samtools merge -@ ${params.threads} ${mergeid}_merged.bam ${bam}
		samtools index -@ ${params.threads} ${mergeid}_merged.bam
		echo 'samtools merge -@ ${params.threads} ${mergeid}_merged.bam ${bam}' > samtools_merge_parameters_${mergeid}.txt
		"""
	}

	singleBams
	.mix(mergedBams)
	.map { mergeid, id, bam, controlid, mark, bam_index ->
		[ mergeid, id, bam, controlid, mark, bam_index ].flatten()
	}
	.set { bamcoverage_mergedbams }

	if (params.mode != "rna" && params.mode != "gro") {
	process bamCoverage {

		publishDir "${params.outdir}/${params.mode}/tracks", mode: 'copy'

		input:
		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bamcoverage_mergedbams

		output:
		file("${mergeid}.RPKMnorm.bw") into bigwigs
		file("bamcoverage_parameters_${mergeid}.txt")
		file("${mergeid}.bamcoverage_report.txt")

		script:
		if (params.mode == "chip" || params.mode == "dnase")
		"""
		bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --smoothLength ${params.bamcoverage_smooth} -e ${params.bamcoverage_e} --centerReads 2> ${mergeid}.bamcoverage_report.txt
		echo 'bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --smoothLength ${params.bamcoverage_smooth} -e ${params.bamcoverage_e} --centerReads' > bamcoverage_parameters_${mergeid}.txt
		"""

		else if (params.mode == "mnase" && params.lib == "s")
		"""
		bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --smoothLength ${params.bamcoverage_smooth} -e ${params.bamcoverage_e} --centerReads --minFragmentLength 130 --maxFragmentLength 200
		"""

		else if (params.mode == "mnase" && params.lib == "p")
		"""
		bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --smoothLength ${params.bamcoverage_smooth} -e ${params.bamcoverage_e} --centerReads --MNase
		"""
		}
	} else {
	process bamCoverage {

		publishDir "${params.outdir}/${params.mode}/tracks", mode: 'copy'

		input:
		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bamcoverage_mergedbams

		output:
		file("${mergeid}.RPKMnorm.fwd.bw") into bigwigs
		file("${mergeid}.RPKMnorm.rev.bw") into bigwigs
		file("bamcoverage_parameters_${mergeid}.txt")
		file("${mergeid}.bamcoverage_report_fwd.txt")
		file("${mergeid}.bamcoverage_report_rev.txt")

		script:
		"""
		bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.fwd.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --filterRNAstrand forward 2> ${mergeid}.bamcoverage_report_fwd.txt
		bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.rev.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --filterRNAstrand reverse 2> ${mergeid}.bamcoverage_report_revtxt
		echo 'bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.fwd.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --filterRNAstrand forward' > bamcoverage_parameters_${mergeid}.txt
		echo 'bamCoverage -b ${bam} -o ${mergeid}.RPKMnorm.rev.bw -of bigwig -bs ${params.bamcoverage_bs} -p ${params.threads} --normalizeUsingRPKM --filterRNAstrand reverse' >> bamcoverage_parameters_${mergeid}.txt
		"""
		}
	}
}

// multiqc general results with bbduk
if (params.multiqc == true && params.bbduk == true && params.lib == "s") {

	process multiqc {

 		publishDir "${params.outdir}/${params.mode}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_multiqc_s.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_multiqc_s.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
	}
}

if (params.multiqc == true && params.bbduk == true && params.lib == "p") {

	process multiqc {

 		publishDir "${params.outdir}/${params.mode}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_multiqc_p.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_multiqc_p.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
	}
}

// multiqc general without bbduk
if (params.multiqc == true && params.bbduk == false && params.lib == "s") {

	process multiqc {

 		publishDir "${params.outdir}/${params.mode}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from pre_fastqc_multiqc_s.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
	}
}

if (params.multiqc == true && params.bbduk == false && params.lib == "p") {

	process multiqc {

 		publishDir "${params.outdir}/${params.mode}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from pre_fastqc_multiqc_p.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
	}
}

///////////////////////////////////////////// DOWNSTREAM ANALYSIS /////////////////////////////////////////////

// chip dnase atac downstream analysis
if (params.downstream_analysis == true && (params.mode == "chip" || params.mode == "dnase" || params.mode == "atac")) {

	treat = Channel.create()
 	control = Channel.create()
 	bam_grouping.choice(treat, control) {
 		it[4] == 'input' ? 1 : 0
 	}

 	process estimate_fragment_size {

 		publishDir "${params.outdir}/${params.mode}/${id}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from spp_bams

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads}
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	// get bams with no control and put them into bams_no_input channel
 	bams.tap{ allBams }
 	.filter{
 		it[3] == '-'
 	}.map {
 		[it[0], it[1], it[2], it[4], it[5], it[6]]
 	}.tap{ bams_no_input }

 	// pair bam files with their control files
 	bams_with_input = control.filter{
 		it[3] != '-'
 	}
 	.cross(allBams) { it[3] }.map{ c, t ->
 		[t[0], t[1], t[2], c[2], t[4], t[5], t[6] ]
 	}

 	// create more bam files with input and no input for downstream analysis
 	bams_with_input.into {
 		macs_input_narrow_bams_noegs
 		macs_input_narrow_bams_egs
 		epic_input_broad_bams_noegs
 		epic_input_broad_bams_egs
 	}

 	bams_no_input.into {
 		macs_no_input_narrow_bams_egs
 		macs_no_input_narrow_bams_noegs
 		macs_no_input_broad_bams
 	}


 	// call peaks for files with input using macs with automatically calculated egs using epic effective
 	if (!params.macs_g && !params.epic_egs) {
 	process narrow_peak_calling_WI {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_WI.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from macs_input_narrow_bams_noegs
 		val egs_size from egs_size_macs_input

 		output:
 		set mergeid, id, file("${id}_peaks.narrowPeak"), mark, fragLen into narrow_peaks_anno_WI
 		file("${id}.macs2_report.txt")
 		file("macs2_parameters_${id}.txt")
 		file("*")

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		if(params.lib == "s" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.lib == "p" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "dnase" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "atac" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad' > macs2_parameters_${id}.txt
 		"""
 	}
 }

 // macs peak calling using manual egs
 	if (params.macs_g == "*" && params.epic_egs == "*") {
 	process narrow_peak_calling_WI {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_WI.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from macs_input_narrow_bams_egs

 		output:
 		set mergeid, id, file("${id}_peaks.narrowPeak"), mark, fragLen into narrow_peaks_anno_WI_m
 		file("${id}.macs2_report.txt")
 		file("macs2_parameters_${id}.txt")
 		file("*")

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		if(params.lib == "s" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.lib == "p" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "dnase" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "atac" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad' > macs2_parameters_${id}.txt
 		"""
 	}
 }

 	// call peaks for files without input using macs
 	if (!params.macs_g && !params.epic_egs) {
 	process narrow_peak_calling_NI {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from macs_no_input_narrow_bams_noegs
 		val egs_size from egs_size_macs_no_input

 		output:
 		set mergeid, id, file("${id}_peaks.narrowPeak"), mark, fragLen into narrow_peaks_anno_NI
 		file("${id}.macs2_report.txt")
 		file("macs2_parameters_${id}.txt")
 		file("*")

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		if(params.lib == "s" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.lib == "p" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAMPE -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${egs_size} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "dnase" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mdoe == "atac" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${egs_size} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad' > macs2_parameters_${id}.txt
 		"""
 	}
 }

 	// call peaks no input with manual egs
 	if (params.macs_g == "*" && params.epic_egs == "*") {
 	process narrow_peak_calling_NI {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from macs_no_input_narrow_bams_egs

 		output:
 		set mergeid, id, file("${id}_peaks.narrowPeak"), mark, fragLen into narrow_peaks_anno_NI_m
 		file("${id}.macs2_report.txt")
 		file("macs2_parameters_${id}.txt")
 		file("*")

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		if(params.lib == "s" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.lib == "p" && params.mode == "chip")
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAMPE -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir . -f BAMPE -g ${params.macs_g} -q ${params.macs_qvalue} -B --SPMR --nomodel --extsize=${fragLen}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mode == "dnase" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize}' > macs2_parameters_${id}.txt
 		"""

 		else if(params.mdoe == "atac" && (params.lib == "s" || params.lib == "p"))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad 2> ${id}.macs2_report.txt
 		echo 'macs2 callpeak -t ${bam} -n ${id} --outdir . -f BAM -g ${params.macs_g} -q ${params.macs_qvalue} --nomodel --extsize 73 --shift 37 --broad' > macs2_parameters_${id}.txt
 		"""
 	}
 }

 	// call broad peaks for files with input using epic automatic egs
	if (params.mode == "chip" && !params.macs_g && !params.epic_egs) {
 	process broad_peak_calling {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_epic.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from epic_input_broad_bams_noegs
 		val egs_ratio from egs_ratio_epic

 		output:
 		set mergeid, id, file("${id}_epic.bed"), mark, fragLen into broad_peaks_anno
 		file("${id}.epic_report.txt")
 		file("epic_parameters_${id}.txt")
 		file("*")

 		script:
 		if (params.lib == "s" && !params.macs_g && !params.epic_egs)
 		"""
 		bedtools bamtobed -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} > ${id}_epic.bed 2> ${id}.epic_report.txt
 		echo 'epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes}' > epic_parameters_${id}.txt
 		"""

 		else if (params.lib == "p" && !params.macs_g && !params.epic_egs)
 		"""
 		bedtools bamtobed -bedpe -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -bedpe -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} --pair-end > ${id}_epic.bed 2> ${id}.epic_report.txt
 		echo 'epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} --pair-end' > epic_parameters_${id}.txt
 		"""
		}
	}

	// epic peak calling using manual egs
	if (params.mode == "chip" && params.macs_g == "*" && params.epic_egs == "*") {
 	process broad_peak_calling {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_epic.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from epic_input_broad_bams_egs

 		output:
 		set mergeid, id, file("${id}_epic.bed"), mark, fragLen into broad_peaks_anno_m
 		file("${id}.epic_report.txt")
 		file("epic_parameters_${id}.txt")
 		file("*")

 		script:
 		if (params.lib == "s" && !params.macs_g && !params.epic_egs)
 		"""
 		bedtools bamtobed -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${params.epic_egs} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} > ${id}_epic.bed 2> ${id}.epic_report.txt
 		echo 'epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${params.epic_egs} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes}' > epic_parameters_${id}.txt
 		"""

 		else if (params.lib == "p" && !params.macs_g && !params.epic_egs)
 		"""
 		bedtools bamtobed -bedpe -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -bedpe -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${params.epic_egs} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} --pair-end > ${id}_epic.bed 2> ${id}.epic_report.txt
 		echo 'epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${params.epic_egs} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.epic_qvalue} -cs ${chromSizes} --pair-end' > epic_parameters_${id}.txt
 		"""
		}
	}
} // closing chip dnase atac

// mnase downstream analysis
if (params.mode == "mnase" && params.downstream_analysis == true) {

	process nucleosome_calling {

 		publishDir "${params.outdir}/${params.mode}/${id}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_danpos.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bam_danpos

 		output:
 		set mergeid, id, file("${id}/pooled/${id}.sorted.mapped.bam.positions.txt"), mark, fragLen into nucleosome_files
 		file("${id}.danpos_report.txt")
 		file("danpos_parameters_${id}.txt")
 		file("*")

 		script:
 		if (params.lib == "s")
 		"""
 		bam2bed --do-not-sort < ${bam} > ${id}.sorted.mapped.bam.bed
 		python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 0 -o ${id} -jw ${params.danpos_jw} -jd ${params.danpos_jd} 2> ${id}.danpos_report.txt
 		mv ${id}/pooled/${id}.sorted.mapped.bam.positions.xls ${id}/pooled/${id}.sorted.mapped.bam.positions.txt
 		echo 'python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 0 -o ${id} -jw ${params.danpos_jw} -jd ${params.danpos_jd}' > danpos_parameters_${id}.txt
 		"""

 		else if (params.lib == "p")
 		"""
 		bam2bed --do-not-sort < ${bam} > ${id}.sorted.mapped.bam.bed
 		python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 1 -o ${id} -jw ${params.danpos_jw} -jd ${params.danpos_jd} 2> ${id}.danpos_report.txt
 		mv ${id}/pooled/${id}.sorted.mapped.bam.positions.xls ${id}/pooled/${id}.sorted.mapped.bam.positions.txt
 		echo 'python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 1 -o ${id} -jw ${params.danpos_jw} -jd ${params.danpos_jd}' > danpos_parameters_${id}.txt
 		"""
 		}
	}
//closing mnase

// rna gro downstream analysis
if (params.downstream_analysis == true && (params.mode == "rna" || params.mode == "gro")) {

	if (!params.strand_info) {
 	exit 1, "Please specify strand information. Available: unstranded, frFirstStrand, frSecondStrand. If you are unsure, run the pipeline using --strandInfo unstranded and --subsample and then look in your qc folder for information on the strandedness of your dataset."
 	}

 	if (!params.dge_file) {
 	exit 1, "Please specify a experiment config file. Check README for exact information. Includes condition information required for differential gene expression analysis."
 	}

 	// set up experiment file
 	exp_file = file(params.dge_file)

 	// quality control with qorts
 	if (params.mode == "rna") {

 		process qorts {

 		publishDir "${params.outdir}/${params.mode}/${id}/qc", mode: 'copy'

 		input:
 		file gtf_file
 		file fasta_file
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_qorts

 		output:
 		file 'QC/*' into qc_files
 		file("${id}.qorts_report.txt")
 		file("qorts_parameters_${id}.txt")

 		script:
 		if (params.lib == "s")
 		"""
 		java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id} --singleEnded ${bam} ${gtf_file} QC 2> ${id}.qorts_report.txt
 		echo 'java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id} --singleEnded ${bam} ${gtf_file} QC' > qorts_parameters_${id}.txt
 		"""

 		if (params.lib == "p")
 		"""
 		java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id} ${bam} ${gtf_file} QC 2> ${id}.qorts_report.txt
 		echo 'java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id} ${bam} ${gtf_file} QC' > qorts_parameters_${id}.txt
 		"""
 	}
 }

 	// preseq for rna data
 	if (params.mode == "rna") {

 		process preseq {

 		publishDir "${params.outdir}/${params.mode}/${id}/preseq", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_preseq

 		output:
 		file("${id}.c_curve.txt") into preseq_results
 		file("${id}.preseq_report.txt")
 		file("preseq_parameters_${id}.txt")

 		script:
 		if (params.lib == "s")
 		"""
 		preseq c_curve -l 9000000000 -B -o ${id}.c_curve.txt ${bam} 2> ${id}.preseq_report.txt
 		echo 'preseq c_curve -l 9000000000 -B -o ${id}.c_curve.txt ${bam}' > preseq_parameters_${id}.txt
 		"""

 		if (params.lib == "p")
 		"""
 		preseq c_curve -l 9000000000 -P -B -o ${id}.c_curve.txt ${bam} 2> ${id}.preseq_report.txt
 		echo 'preseq c_curve -l 9000000000 -P -B -o ${id}.c_curve.txt ${bam}' > preseq_parameters_${id}.txt
 		"""
 		}
 	}

 	// read counting using featurecounts
 	process read_counting {

 		publishDir "${params.outdir}/${params.mode}/${id}/featurecounts", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_featurecounts
 		file gtf_file

 		output:
 		file("${id}_gene.featureCounts.txt") into geneCounts
 		file("${id}_gene.featureCounts.txt.summary") into featureCounts_logs
 		file("${id}_biotype_counts.txt") into featureCounts_biotype
 		file("${id}.featurecounts_report.txt")
 		file("featurecounts_parameters_${id}.txt")

 		script:
 		if (params.strand_info == "unstranded" && params.lib == "s")
 		"""
 		featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 0 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam}' > featurecounts_parameters_${id}.txt
 		"""

 		else if (params.strand_info == "frFirstStrand" && params.lib == "s")
 		"""
 		featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 2 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam}' > featurecounts_parameters_${id}.txt
 		"""

 		else if (params.strand_info == "frSecondStrand" && params.lib == "s")
 		"""
 		featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 1 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam}' > featurecounts_parameters_${id}.txt
 		"""

 		else if (params.strand_info == "unstranded" && params.lib == "p")
 		"""
 		featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 0 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam}' > featurecounts_parameters_${id}.txt
 		"""

 		else if (params.strand_info == "frFirstStrand" && params.lib == "p")
 		"""
 		featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 2 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam}' > featurecounts_parameters_${id}.txt
 		"""

 		else if (params.strand_info == "frSecondStrand" && params.lib == "p")
 		"""
 		featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam} 2> ${id}.featurecounts_report.txt
 		featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 1 ${bam}
 		cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 		echo 'featureCounts -a ${gtf_file} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam}' > featurecounts_parameters_${id}.txt
 		"""
 	}

 	// rpkm and tpm scores using stringtie
 	process stringtie {

 		publishDir "${params.outdir}/${params.mode}/${id}/stringtie", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_stringtie
 		file gtf_file

 		output:
 		file("${id}.gene_abundance.txt")
 		file("${id}.cov_refs.gtf")
 		file("stringtie_log") into stringtie_log

 		script:
 		"""
 		stringtie ${bam} -G ${gtf_file} -A ${id}.gene_abundance.txt -C ${id}.cov_refs.gtf -e -b ${id}_ballgown -p ${params.threads} > stringtie_log
 		"""
 	}

 	// multiqc
 	if (params.multiqc == true) {

 		process multiqc {

 		publishDir "${params.outdir}/${params.mode}/${id}/multiqc", mode: 'copy'

 		input:
 		file ('preseq/*') from preseq_results.flatten().toList()
 		file ('featurecounts/*') from featureCounts_logs.flatten().toList()
 		file ('stringtie/*') from stringtie_log.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 		}
 	}

 	// dge
 	if (params.dge == true) {

 		process dge {

 		publishDir "${params.outdir}/${params.mode}/dge", mode: 'copy'

 		input:
 		file input_files from geneCounts.toSortedList()
 		file exp_file

 		output:
 		file "*" into dge_results

 		script:
 		"""
 		Rscript $baseDir/bin/dge.R ${exp_file} $input_files
 		"""
 		}
 	}
} // closing rna gro

// analysis mode
if (params.mode == "analysis") {

	if (params.function == "predictEnhancers") {

		// Parse config file
 		bedgraphs = Channel
 		.from(config_file.readLines())
 		.map { line ->
 			list = line.split()
 			id = list[0]
 			bedgraph = file(list[1])
 			[ id, bedgraph ]
 		}

 		// Step 1. FASTA index
 		process fasta_index {

 			input:
 			file fasta_file

 			output:
 			file "chromSizes.txt" into chromSizes

 			script:
 			"""
 			samtools faidx ${fasta_file}
			awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
			"""
 		}

 		// Step 2. Split genome into 200bp windows
 		process generate_genome_bins {

 			input:
 			file chromSizes from chomSizes.val

 			output:
 			file("split_genome.txt") into genome

 			script:
 			"""
 			bedtools makewindows -w 200 -g ${chromSizes} | sort-bed - | awk -v OFS='\t' '{if (\$1 ~ /chr/) print \$0}' > split_genome.txt
 			"""
 		}

 		genome.into {
 			genome1
 			genome2
 		}

 		// Step 3. Convert bedgraphs into beds
 		process bedgraph_to_bed {

			input:
			set id, file(bedgraph) from bedgraphs

			output:
			set id, file("${id}.bed") into beds

			script:
			"""
			awk '{print \$1"\t"\$2"\t"\$3"\t.\t"\$4}' ${bedgraph} | sort-bed - | awk -v OFS='\t' '{if (\$1 ~ /chr/) print \$0}' > ${id}.bed
			"""
		}

		// Step 4. Calculate Mean Coverage
		process mean_coverage_over_genome {

			input:
			set id, file(bed) from beds
			file genome from genome1.val

			output:
			set id, file("${id}.coverage_only.txt") into coverage_files

			script:
			"""
			bedmap --echo --mean --delim '\t' ${genome} ${bed} > ${id}.coverages.txt
			awk -v OFS='\t' '{print \$4}' ${id}.coverages.txt > ${id}.coverage_only.txt
			"""
		}

		coverage_files.into {
		coverage_files_1
		coverage_files_2
		coverage_files_3
		coverage_files_4
		}

		dnase = Channel.create()
		dnase = coverage_files_1.filter {
			it[0] == "DNase"
		}

		h3k1 = Channel.create()
		h3k1 = coverage_files_2.filter {
			it[0] == "H3K4me1"
		}

		h3k27ac = Channel.create()
		h3k27ac = coverage_files_3.filter {
			it[0] == "H3K27Ac"
		}

		h3k3 = Channel.create()
		h3k3 = coverage_files_4.filter {
			it[0] == "H3K4me3"
		}

		// Step 5. Combine coverages into one final dataset
		process merge_coverages {

			input:
			set id, file(dnase_file) from dnase
			set id, file(h3k1_file) from h3k1
			set id, file(h3k27ac_file) from h3k27ac
			set id, file(h3k3_file) from h3k3
			file genome from genome2.val

			output:
			file("final.dataset.txt") into final_dataset

			script:
			"""
			paste ${genome} ${dnase_file} ${h3k27ac_file} ${h3k1_file} ${h3k3_file} > final.dataset.txt
			"""
		}

		// Step 6: Generate list of putative enhancers
		process predict_enhancers {

			input:
			file(data) from final_dataset

			output:
			file("predicted.enhancers.txt") into enhancers

			script:
			"""
			Rscript '$baseDir/bin/enhancer_prediction.R' ${data} '$baseDir/bin/mlModel.rds'
			"""
		}

		// Step 7. Merge putative enhancers
		process stitching {

			publishDir "${params.outdir}/analysis", mode: 'copy'

			input:
			file(enhancer_file) from enhancers

			output:
			file("predicted_enhancers.txt")

			script:
			"""
			sort-bed ${enhancer_file} | bedops --merge - | awk -v OFS='\t' '{print \$1, \$2, \$3, "Enhancer.Prediction." NR}' > predicted_enhancers.txt
			"""
		}
	} // closing predictEnhancers

	if (params.function == "geneExpressionNearPeaks") {

		feature_file = file(params.feature_file)
 		stringtie_file = file(params.stringtie_file)

 		process calculate_gene_expression_near_peaks {

 			publishDir "${params.outdir}/analysis", mode: 'copy'

 			input:
 			gtf_file
 			peak_file
 			stringtie_file

 			output:
 			file("geneExpression_near_ChIPseq_peaks.txt")

 			script:
 			"""
 			Rscript '$baseDir/bin/calculate_gene_expression_near_peaks.R' ${gtf_file} ${feature_file} ${stringtie_file}
 			"""
 		}
	} // closing geneExpressionNearPeaks
} // closing analysis

// on completion
 workflow.onComplete {
 	println ""
 	println "Workflow completed on: $workflow.complete"
 	println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
 	println "Workflow Duration: $workflow.duration"
 	println ""
 	println "Submit issues/requests on GitHub or e-mail cag104@ucsd.edu"
 	println ""
 }