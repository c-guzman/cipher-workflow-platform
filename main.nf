#!/usr/bin/env nextflow
// Author: Carlos Guzman
// E-mail: cag104@ucsd.edu

/*
 * DEFINE DEFAULT PARAMETERS
 */

 params.threads = 1
 params.minid = 0.76
 params.help = false
 params.subsample = false
 params.outdir = './results'
 params.qvalue = 0.05
 params.epic_w = 200
 params.epic_g = 3
 params.maxindel = '200k'
 params.intronlen = 20
 params.egs = false
 params.egs_ratio = false
 params.aligner = 'bbmap'
 //params.binSize = 10
 //params.smoothLen = 50
 params.aligner = 'bbmap'

 // bwa defaults
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

 // bowtie2 defaults
 params.bt2_D = 20
 params.bt2_R = 3
 params.bt2_N = 0
 params.bt2_L = 20
 params.bt2_i = '5,1,0.50'
 params.bt2_trim5 = 0
 params.bt2_trim3 = 0
 params.local = false

 // hisat2 defaults
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

 // star defaults
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

 // PRINT HELP
 if (params.help) {
 	log.info ''
 	log.info '~ C I P H E R ~ Version 1.1.0'
 	log.info '********************************************************************************************************************************************************'
 	log.info ''
 	log.info 'REQUIRED PARAMETERS:'
 	log.info '===================='
 	log.info '--mode				Choose from available: chip, rna, gro, mnase, dnase, atac, analysis.'
 	log.info '--config			Configuration file with sample information. Check README for more information.'
 	log.info '--fasta				Reference genome in FASTA format.'
 	log.info '--gtf				Reference genome in GTF format.'
 	log.info '--lib				Library information. "s" for single-stranded data, "p" for pair-ended data.'
 	log.info '--readLen			The length of your reads.'
 	log.info ''
 	log.info 'RNA-seq MODE ONLY:'
 	log.info '=================='
 	log.info '--strandInfo			Strandedness information. Choose from "unstranded", "frFirstStrand", or "frSecondStrand".'
 	log.info '--expInfo			Experiment config file for RNA-seq data DGE analysis. Check README for more information.'
 	log.info ''
 	log.info 'ANALYSIS MODE ONLY:'
 	log.info '==================='
 	log.info '--analysis			Choose from available: "predictEnhancers".'
 	log.info ''
 	log.info 'OPTIONAL PARAMETERS:'
 	log.info '===================='
 	log.info '--threads			Number of threads. (Default: 1)'
 	log.info '--aligner			The alignment software your workflow will use. Available: bbmap, bowtie2, bwa, hisat2, star. (Default: bbmap)'
 	log.info '--minid				Minimum alignment identity to look for during BBMap mapping. Higher is faster and less sensitive. (Default: 0.76)'
 	log.info '--qvalue			Minimum FDR cutoff for peak detection in MACS2 and EPIC. (Default: 0.05)'
 	log.info '--epic_w			Size of the windows used to scan the genome for peak detection in EPIC. (Default: 200)'
 	log.info '--epic_g			A multiple of epic_w used to determine the gap size in EPIC. (Default: 3)'
 	log.info '--maxindel			Maximum indel length searched during mapping. 200k recommended for vertebrate genomes. (Default: 200k)'
 	log.info '--intronlen			Maximum intron length during mapping. 20 recommended for vertebrate genomes. (Default: 20)'
 	log.info '--egs				The effective genome size of your species. (Default: Automatically calculated - requires approximately 80GB of RAM)'
 	log.info '--egs_ratio			Effective genome as fraction of the genome size. Must be between 0 and 1. Check EPIC GitHub for more information. (Default: Automatically calculated - requires approximately 80GB of RAM)'
 	log.info '--outdir			Name of output directory. (Default: results)'
 	log.info ''
 	log.info 'ALIGNER-SPECIFIC OPTIONAL PARAMETERS:'
 	log.info '====================================='
 	log.info ''
 	log.info 'BWA:'
 	log.info '--bwa_k				See -k option in BWA user manual for more information. (Default: 19)'
 	log.info '--bwa_w				See -w option in BWA user manual for more information. (Default: 100)'
 	log.info '--bwa_d				See -d option in BWA user manual for more information. (Default: 100)'
 	log.info '--bwa_r				See -r option in BWA user manual for more information. (Default: 1.5)'
 	log.info '--bwa_c				See -c option in BWA user manual for more information. (Default: 10000)'
 	log.info '--bwa_A				See -A option in BWA user manual for more information. (Default: 1)'
 	log.info '--bwa_B				See -B option in BWA user manual for more information. (Default: 4)'
 	log.info '--bwa_O				See -O option in BWA user manual for more information. (Default: 6)'
 	log.info '--bwa_E				See -E option in BWA user manual for more information. (Default: 1)'
 	log.info '--bwa_L				See -L option in BWA user manual for more information. (Default: 5)'
 	log.info '--bwa_U				See -U option in BWA user manual for more information. (Default: 9)'
 	log.info '--bwa_T				See -T option in BWA user manual for more information. (Default: 30)'
 	log.info ''
 	log.info 'BOWTIE2:'
 	log.info '--bt2_D				See -D option in Bowtie2 user manual for more information. (Default: 20)'
 	log.info '--bt2_R				See -R option in Bowtie2 user manual for more information. (Default: 3)'
 	log.info '--bt2_N				See -N option in Bowtie2 user manual for more information. (Default: 0)'
 	log.info '--bt2_L				See -L option in Bowtie2 user manual for more information. (Default: 20)'
 	log.info '--bt2_i				See -i option in Bowtie2 user manual for more information. (Default: 5,1,0.50)'
 	log.info '--bt2_trim5			See --trim5 option in Bowtie2 user manual for more information. (Default: 0)'
 	log.info '--bt2_trim3			See --trim3 option in Bowtie2 user manual for more information. (Default: 0)'
 	log.info '--local				Set this parameter to map alignments in local mode. (Default: false)'
 	log.info ''
 	log.info 'HISAT2:'
 	log.info '--hs_k				See -k option in HISAT2 user manual for more information. (Default: 5)'
 	log.info '--hs_trim5			See --trim5 option in HISAT2 user manual for more information. (Default: 0)'
 	log.info '--hs_trim3			See --trim3 option in HISAT2 user manual for more information. (Default: 0)'
 	log.info '--hs_mp				See --mp option in HISAT2 user manual for more information. (Default: 6,2)'
 	log.info '--hs_sp				See --sp option in HISAT2 user manual for more information. (Default: 2,1)'
 	log.info '--hs_np				See --np option in HISAT2 user manual for more information. (Default: 1)'
 	log.info '--hs_rdg			See --rdg option in HISAT2 user manual for more information. (Default: 5,3)'
 	log.info '--hs_rfg			See --rfg option in HISAT2 user manual for more information. (Default: 5,3)'
 	log.info '--hs_pen_cansplice		See --pen-cansplice option in HISAT2 user manual for more information. (Default: 0)'
 	log.info '--hs_pen_noncansplice		See --pen-nonansplice option in HISAT2 user manual for more information. (Default: 12)'
 	log.info '--hs_min_intronlen		See --min-intronlen option in HISAT2 user manual for more information. (Default: 20)'
 	log.info '--hs_max_intronlen		See --max-intronlen option in HISAT2 user manual for more information. (Default: 500000)'
 	log.info '--hs_max_seeds			See --max-seeds option in HISAT2 user manual for more information. (Default: 5)'
 	log.info ''
 	log.info 'STAR:'
 	log.info '--star_clip3pNbases		See --clip3pNbases option in STAR user manual for more information. (Default: 0)'
 	log.info '--star_clip5pNbases		See --clip5pNbases option in STAR user manual for more information. (Default: 0)'
 	log.info '--star_outFilterMultimapScoreRange	See --outFilterMultimapScoreRange option in STAR user manual for more information. (Default: 1)'
 	log.info '--star_outFilterMultimapNmax		See --outFilterMultimapNmax option in STAR user manual for more information. (Default: 10)'
 	log.info '--star_outFilterMismatchNmax		See --outFilterMismatchNmax option in STAR user manual for more information. (Default: 10)'
 	log.info '--star_outFilterScoreMin		See --outFilterScoreMin option in STAR user manual for more information. (Default: 0)'
 	log.info '--star_alignEndsType			See --alignEndsType option in STAR user manual for more information. (Default: Local)'
 	log.info '--star_winAnchorMultimapNmax		See --winAnchorMultimapNmax option in STAR user manual for more information. (Default: 50)'
 	log.info '--star_quantMode			See --quantMode option in STAR user manual for more information. (Default: -)'
 	log.info '--star_twopassMode			See --twopassMode option in STAR user manual for more information. (Default: None)'
 	log.info ''
 	log.info 'FOR TESTING AND OPTIMIZING:'
 	log.info '==========================='
 	log.info '--subsample			Set this flag to subsample reads for testing.'
 	log.info ''
 	log.info '********************************************************************************************************************************************************'
 	exit 1
 }

 // PARAMETER CHECKS
 if (!params.config) {
 	exit 1, "Please specify a config file"
 }

 if (!params.fasta) {
 	exit 1, "Please specify a FASTA file"
 }

 if (!params.gtf) {
 	exit 1, "Please specify a GTF file"
 }

 if (!params.readLen) {
 	exit 1, "Please specify the length of your sequenced reads before trimming"
 }

 if (!params.mode) {
 	exit 1, "Please specify a pipeline mode. Available: chip, rna, mnase, dnase, analysis"
 }

 if (!params.lib) {
 	exit 1, "Please specify the strandedness of your experiment. 's' for single-ended and 'p' for pair-ended data"
 }

 if (params.mode != 'analysis' && params.analysis) {
 	exit 1, "--analysis should only be used when running analysis mode (--mode analysis)"
 }

 if (params.aligner != 'bbmap' || params.aligner != 'bowtie2' || params.aligner != 'bwa' || params.aligner != 'hisat2' || params.aligner != 'star') {
 	exit 1, "--aligner should be one of the following: bbmap, bowtie2, bwa, hisat2, star"
 }

 // SET UP FILES
 config_file = file(params.config)
 fasta_file = file(params.fasta)
 gtf_file = file(params.gtf)

 // DISPLAY RUN INFORMATION
 log.info ''
 log.info 'VERSION 1.0.0 - CIPHER'
 log.info '**********************************************'
 log.info ''
 log.info "mode:		${params.mode}"
 log.info "config:		${config_file}"
 log.info "fasta:		${fasta_file}"
 log.info "gtf:		${gtf_file}"
 log.info "read-length:	${params.readLen}"
 log.info "library:	${params.lib}"
 log.info ""
 log.info "aligner:		${params.aligner}"
 log.info ''
 log.info "subsample:	${params.subsample}"
 log.info ''
 log.info '**********************************************'
 log.info ''

 // Download Data

 // SE CHIP
 if (params.mode == 'chip' && params.lib == 's') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	if (!params.egs && !params.egs_ratio) {
 		process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file
 		file "egs_ratio.txt" into egs_ratio

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_WI
 			egs_size_deeptools_NI
 			egs_size_macs_WI
 			egs_size_macs_NI
 		}

 		egs_ratio.map{ file ->
 			file.text.trim() } .set {
 				egs_ratio
 			}

 			egs_ratio.into {
 				egs_ratio_epic_WI
 			}
 		}


 	if (params.aligner == 'bbmap') {
 	   // Generate BBMap Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 if (params.aligner == 'bowtie2') {
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

if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	treat = Channel.create()
 	control = Channel.create()
 	bam_grouping.choice(treat, control) {
 		it[4] == 'input' ? 1 : 0
 	}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_spp

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	// GET BAMS WITH NO CONTROL AND PUT THEM INTO BAMSNOINPUT CHANNEL
 	bams.tap{ allBams }
 	.filter{
 		it[3] == '-'
 	}.map {
 		[it[0], it[1], it[2], it[4], it[5], it[6]]
 	}.tap{ bams_no_input }

 	// PAIR BAMS WITH CONTROLS
 	bams_with_input = control.filter{
 		it[3] != '-'
 	}
 	.cross(allBams) { it[3] }.map{ c, t ->
 		[t[0], t[1], t[2], c[2], t[4], t[5], t[6]]
 	}

 	// CREATE MORE BAM CHANELS WITH INPUT
 	bams_with_input.into {
 		bamsWI_bigwigs
 		bamsWI_macs
 		bamsWI_epic
 		bamsWI_chipqc
 	}

 	// CREATE MORE BAM CHANNELS WITH NO INPUT
 	bams_no_input.into {
 		bamsNI_bigwigs
 		bamsNI_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH INPUT
 	process create_coverage_tracks_WI {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_bigwigs
 		val egs_size from egs_size_deeptools_WI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_WI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks_NI {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bamsNI_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_NI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 WITH INPUT
 	process call_binding_sites_WI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_WI.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_macs
 		val egs_size from egs_size_macs_WI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_WI_anno, np_WI_motifs, np_WI_chipqc

 		script:
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir ${id} -f BAM -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111
 		"""

 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites_NI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bamsNI_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno
 		set mergeid, id, file("${id}/${id}_peaks.broadPeak"), mark, fragLen into bp_NI_anno

 		script:
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAM -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAM -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111 --broad --broad-cutoff ${params.qvalue}
 		"""
 	}

 	// STEP 8 BROAD PEAK CALLING WITH EPIC WITH INPUT
 	process call_broad_binding_sites_WI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_epic.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_epic
 		val egs_ratio from egs_ratio_epic_WI

 		output:
 		set mergeid, id, file("${id}_epic.broadPeak"), mark, fragLen into bp_WI_anno

 		script:
 		"""
 		bedtools bamtobed -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.qvalue} -cs ${chromSizes} > ${id}_epic.broadPeak
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_WI_anno
 		set mergeid, id, file(bp), mark, fragLen from bp_WI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")
 		file("${id}_broadPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${bp} ${fasta_file} -gtf ${gtf_file} > ${id}_broadPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 IDENTIFY POTENTIAL MOTIFS
 	process identify_binding_motifs {

 		publishDir "${params.outdir}/motifs", mode: 'copy'

 		input:
 		set mergeid, id, file(np), mark, fragLen from np_WI_motifs
 		file fasta_file

 		output:
 		file("${id}_motifs/*")

 		script:
 		"""
 		sort -nk8 ${np} | awk '{ a[NR] = \$0 } END { for (i = 1; i <= NR / 10; ++i) print a[i] }' > ${id}.sorted.peak.files.txt
 		perl $baseDir/bin/homer/bin/findMotifsGenome.pl ${id}.sorted.peak.files.txt ${fasta_file} ${id}_motifs -size 100
 		"""
 	}

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	// STEP 12 CHIPQC WITH INPUT
 	process chipqc {

 		publishDir "${params.outdir}/chipqc", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_chipqc
 		set mergeid, id, file(np), mark, fragLen from np_WI_chipqc

 		output:
 		file("${id}/*")

 		script:
 		"""
 		Rscript $baseDir/bin/chipqc.R ${bam} ${np} ${id}
 		"""
 	}

 } // closing bracket SE chip

 // PE CHIP
 if (params.mode == 'chip' && params.lib == 'p') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file
 		file "egs_ratio.txt" into egs_ratio

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_WI
 			egs_size_deeptools_NI
 			egs_size_macs_WI
 			egs_size_macs_NI
 		}

 		egs_ratio.map{ file ->
 			file.text.trim() } .set {
 				egs_ratio
 			}

 			egs_ratio.into {
 				egs_ratio_epic_WI
 			}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	treat = Channel.create()
 	control = Channel.create()
 	bam_grouping.choice(treat, control) {
 		it[4] == 'input' ? 1 : 0
 	}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_pp

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	// GET BAMS WITH NO CONTROL AND PUT THEM INTO BAMSNOINPUT CHANNEL
 	bams.tap{ allBams }
 	.filter{
 		it[3] == '-'
 	}.map {
 		[it[0], it[1], it[2], it[4], it[5], it[6]]
 	}.tap{ bams_no_input }

 	// PAIR BAMS WITH CONTROLS
 	bams_with_input = control.filter{
 		it[3] != '-'
 	}
 	.cross(allBams) { it[3] }.map{ c, t ->
 		[t[0], t[1], t[2], c[2], t[4], t[5], t[6], c[6]]
 	}

 	// CREATE MORE BAM CHANELS WITH INPUT
 	bams_with_input.into {
 		bamsWI_bigwigs
 		bamsWI_macs
 		bamsWI_epic
 		bamsWI_chipqc
 	}

 	// CREATE MORE BAM CHANNELS WITH NO INPUT
 	bams_no_input.into {
 		bamsNI_bigwigs
 		bamsNI_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH INPUT
 	process create_coverage_tracks_WI {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_bigwigs
 		val egs_size from egs_size_deeptools_WI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_WI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks_NI {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bamsNI_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_NI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 WITH INPUT
 	process call_binding_sites_WI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_WI.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_macs
 		val egs_size from egs_size_macs_WI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_WI_anno, np_WI_motifs, np_WI_chipqc

 		script:
 		"""
 		macs2 callpeak -t ${bam} -c ${control} -n ${id} --outdir ${id} -f BAMPE -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111
 		"""

 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites_NI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bamsNI_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno
 		set mergeid, id, file("${id}/${id}_peaks.broadPeak"), mark, fragLen into bp_NI_anno

 		script:
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAMPE -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAMPE -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --seed 111 --broad --broad-cutoff ${params.qvalue}
 		"""
 	}

 	// STEP 8 BROAD PEAK CALLING WITH EPIC WITH INPUT
 	process call_broad_binding_sites_WI {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_epic.val
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_epic
 		val egs_ratio from egs_ratio_epic_WI

 		output:
 		set mergeid, id, file("${id}_epic.broadPeak"), mark, fragLen into bp_WI_anno

 		script:
 		"""
 		bedtools bamtobed -bedpe -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -bedpe -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_ratio} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.qvalue} -cs ${chromSizes} --pair-end > ${id}_epic.bed
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_WI_anno
 		set mergeid, id, file(bp), mark, fragLen from bp_WI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")
 		file("${id}_broadPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${bp} ${fasta_file} -gtf ${gtf_file} > ${id}_broadPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 IDENTIFY POTENTIAL MOTIFS
 	process identify_binding_motifs {

 		publishDir "${params.outdir}/motifs", mode: 'copy'

 		input:
 		set mergeid, id, file(np), mark, fragLen from np_WI_motifs
 		file fasta_file

 		output:
 		file("${id}_motifs/*")

 		script:
 		"""
 		sort -nk8 ${np} | awk '{ a[NR] = \$0 } END { for (i = 1; i <= NR / 10; ++i) print a[i] }' > ${id}.sorted.peak.files.txt
 		perl $baseDir/bin/homer/bin/findMotifsGenome.pl ${id}.sorted.peak.files.txt ${fasta_file} ${id}_motifs -size 100
 		"""
 	}

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	// STEP 12 CHIPQC WITH INPUT
 	process chipqc {

 		publishDir "${params.outdir}/chipqc", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index) from bamsWI_chipqc
 		set mergeid, id, file(np), mark, fragLen from np_WI_chipqc

 		output:
 		file("${id}/*")

 		script:
 		"""
 		Rscript $baseDir/bin/chipqc.R ${bam} ${np} ${id}
 		"""
 	}

 } // closing bracket PE chip

 // SE DNASE
 if (params.mode == 'dnase' && params.lib == 's') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno, np_NI_motifs

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAM -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} --seed 111
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 IDENTIFY POTENTIAL MOTIFS
 	process identify_binding_motifs {

 		publishDir "${params.outdir}/motifs", mode: 'copy'

 		input:
 		set mergeid, id, file(np), mark, fragLen from np_NI_motifs
 		file fasta_file

 		output:
 		file("${id}_motifs/*")

 		script:
 		"""
 		sort -nk8 ${np} | awk '{ a[NR] = \$0 } END { for (i = 1; i <= NR / 10; ++i) print a[i] }' > ${id}.sorted.peak.files.txt
 		perl $baseDir/bin/homer/bin/findMotifsGenome.pl ${id}.sorted.peak.files.txt ${fasta_file} ${id}_motifs -size 100
 		"""
 	}

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket SE dnase

 	// PE DNASE
 if (params.mode == 'dnase' && params.lib == 'p') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno, np_NI_motifs

 		script:
 		shiftsize = Math.round(-1 * ((fragLen as int)/2))
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAMPE -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize=${fragLen} --shift=${shiftsize} --seed 111
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 IDENTIFY POTENTIAL MOTIFS
 	process identify_binding_motifs {

 		publishDir "${params.outdir}/motifs", mode: 'copy'

 		input:
 		set mergeid, id, file(np), mark, fragLen from np_NI_motifs
 		file fasta_file

 		output:
 		file("${id}_motifs/*")

 		script:
 		"""
 		sort -nk8 ${np} | awk '{ a[NR] = \$0 } END { for (i = 1; i <= NR / 10; ++i) print a[i] }' > ${id}.sorted.peak.files.txt
 		perl $baseDir/bin/homer/bin/findMotifsGenome.pl ${id}.sorted.peak.files.txt ${fasta_file} ${id}_motifs -size 100
 		"""
 	}

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket PE dnase

 	// SE MNASE
 	if (params.mode == 'mnase' && params.lib == 's') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap')  {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 1 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 20 -e ${fragLen} --ignoreDuplicates --centerReads --minFragmentLength 130 --maxFragmentLength 200
 		"""
 	}

 	// STEP 7 NUCLEOSOME PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs

 		output:
 		set mergeid, id, file("${id}/pooled/${id}.sorted.mapped.bam.positions.txt"), mark, fragLen into np_NI_anno

 		script:
 		"""
 		bam2bed --do-not-sort < ${bam} > ${id}.sorted.mapped.bam.bed
 		python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 0 -o ${id} --frsz ${fragLen} -jw 40 -jd 150
 		mv ${id}/pooled/${id}.sorted.mapped.bam.positions.xls ${id}/pooled/${id}.sorted.mapped.bam.positions.txt
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_nucleosomes.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_nucleosomes.annotated.txt
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket SE mnase

 	// PE MNASE
 	if (params.mode == 'mnase' && params.lib == 'p') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 1 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 20 -e ${fragLen} --ignoreDuplicates --centerReads --MNase
 		"""
 	}

 	// STEP 7 NUCLEOSOME PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs

 		output:
 		set mergeid, id, file("${id}/pooled/${id}.sorted.mapped.bam.positions.txt"), mark, fragLen into np_NI_anno

 		script:
 		"""
 		bam2bed --do-not-sort < ${bam} > ${id}.sorted.mapped.bam.bed
 		python $baseDir/bin/danpos/danpos.py dpos ${id}.sorted.mapped.bam.bed -m 1 -o ${id} --frsz ${fragLen} -jw 40 -jd 150
 		mv ${id}/pooled/${id}.sorted.mapped.bam.positions.xls ${id}/pooled/${id}.sorted.mapped.bam.positions.txt
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_nucleosomes.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_nucleosomes.annotated.txt
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket PE mnase

 	// SE GRO
 	if (params.mode == 'gro' && params.lib == 's') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs_fwd
 		bams_bigwigs_rev
 		bams_merge
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_fwd_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs_fwd
 		val egs_size from egs_size_deeptools_fwd

 		output:
 		file("${id}.RPGCnorm.fwd.bigWig") into fwd_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --filterRNAstrand forward --Offset 1
 		"""
 	}

 	// STEP 7 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_rev_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, fragLen, file(bam_index) from bams_bigwigs_rev
 		val egs_size from egs_size_deeptools_rev

 		output:
 		file("${id}.RPGCnorm.rev.bigWig") into rev_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.rev.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --filterRNAstrand reverse --Offset 1
 		"""
 	}

 	// Merge replicates
 	singleBams = Channel.create()
 	groupedBams = Channel.create()

 	bams_merge.groupTuple(by: [0,3,4])
 	.choice(singleBams, groupedBams) {
 		it[2].size() > 1 ? 1 : 0
 	}

 	// STEP 8 MERGE REPLICATES
 	process merge_replicates {

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_merge

 		output:
 		set mergeid, id, file(""), mark, fragLen into mergedBams

 		script:
 		def id = id.sort().join(':')
 		"""
 		(
 		samtools view -H ${bam} | grep -v '@SQ';
 		for f in ${bam}; do
 		  samtools view -H \$f | grep '@SQ';
 		done
 		) > header.txt && \
 		samtools merge -@ ${params.threads} -h header.txt ${mergeid}_primary.bam ${bam}
 		"""
 	}

 	singleBams
 	.mix(mergedBams)
 	.map { mergeid, id, bam, controlid, mark, fragLen ->
 		[ mergeid, bam, controlid, mark, fragLen ].flatten()
 	}
 	.into { bams_grohmm }

 	// STEP 9 GROHMM WORKFLOW
 	process groHMM {

 		publishDir "${params.outdir}/groHMM", mode: 'copy'

 		input:
 		set mergeid, file(bam), controlid, mark, fragLen from bams_grohmm

 		output:
 		set mergeid, file("${mergeid}_transcripts.txt") into grohmm_transcripts

 		script:
 		"""
 		Rscript $baseDir/bin/grohmm.R ${bam} ${mergeid} ${params.threads}
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket SE gro

 	// PE GRO
 	if (params.mode == 'gro' && params.lib == 'p') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }


	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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
 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs_fwd
 		bams_bigwigs_rev
 		bams_merge
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_fwd_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs_fwd
 		val egs_size from egs_size_deeptools_fwd

 		output:
 		file("${id}.RPGCnorm.fwd.bigWig") into fwd_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --filterRNAstrand forward --Offset 1
 		"""
 	}

 	// STEP 7 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_rev_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, fragLen, file(bam_index) from bams_bigwigs_rev
 		val egs_size from egs_size_deeptools_rev

 		output:
 		file("${id}.RPGCnorm.rev.bigWig") into rev_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.rev.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --filterRNAstrand reverse --Offset 1
 		"""
 	}

 	// Merge replicates
 	singleBams = Channel.create()
 	groupedBams = Channel.create()

 	bams_merge.groupTuple(by: [0,3,4])
 	.choice(singleBams, groupedBams) {
 		it[2].size() > 1 ? 1 : 0
 	}

 	// STEP 8 MERGE REPLICATES
 	process merge_replicates {

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_merge

 		output:
 		set mergeid, id, file(""), mark, fragLen into mergedBams

 		script:
 		def id = id.sort().join(':')
 		"""
 		(
 		samtools view -H ${bam} | grep -v '@SQ';
 		for f in ${bam}; do
 		  samtools view -H \$f | grep '@SQ';
 		done
 		) > header.txt && \
 		samtools merge -@ ${params.threads} -h header.txt ${mergeid}_primary.bam ${bam}
 		"""
 	}

 	singleBams
 	.mix(mergedBams)
 	.map { mergeid, id, bam, controlid, mark, fragLen ->
 		[ mergeid, bam, controlid, mark, fragLen ].flatten()
 	}
 	.into { bams_grohmm }

 	// STEP 9 GROHMM WORKFLOW
 	process groHMM {

 		publishDir "${params.outdir}/groHMM", mode: 'copy'

 		input:
 		set mergeid, file(bam), controlid, mark, fragLen from bams_grohmm

 		output:
 		set mergeid, file("${mergeid}_transcripts.txt") into grohmm_transcripts

 		script:
 		"""
 		Rscript $baseDir/bin/grohmm.R ${bam} ${mergeid} ${params.threads}
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket PE gro

 	// SE RNA
 	if (params.mode == 'rna' && params.lib == 's') {

 		if (!params.strandInfo) {
 	exit 1, "Please specify strand information. Available: unstranded, frFirstStrand, frSecondStrand. If you are unsure, run the pipeline using --strandInfo unstranded and --subsample and then look in your qc folder for information on the strandedness of your dataset."
 }

 		if (!params.expInfo) {
 	exit 1, "Please specify a experiment config file. Check README for exact information. Typically includes condition information for differential gene expression analysis."
 }

 exp_file = file(params.expInfo)

 	// Parse config file
 	fastqs = Channel
 	.from(config_file.readLines())
 	.map { line ->
 			list = line.split()
 			mergeid = list[0]
 			id = list[1]
 			path1 = file(list[2])
 			[ mergeid, id, path1 ]
 		}

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1) from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz") into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	if (!params.egs) {
 		process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_fwd
 			egs_size_deeptools_rev
 		}
 	}

 	// Generate BBMap Index
 	if (params.aligner == 'star') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

  	if (params.aligner == 'hisat2') {
	// Generate HISAT2 Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/hisat2_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("*") into hisat22_index

 		script:
 		"""
 		hisat2-build -p ${params.threads} ${fasta_file} genome
 		"""
 	}
}

	if (params.aligner == 'star') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz") into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs, qorts_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 intronlen=${params.intronlen} maxindel=${params.maxindel} ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 // STEP 4 MAPPING WITH HISAT2
if (params.aligner == 'hisat2') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from bbmap_trimmed_fastqs
 		file("*") from hisat2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond
 		file("${id}.hs2.metricsFile.txt")

 		script:
 		"""
 		hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 if (params.aligner == 'star') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from bbmap_trimmed_fastqs
 		file("indexFiles/*") from star_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond

 		script:
 		"""
 		STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode}
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}Aligned.sortedByCoord.out.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// STEP 5 CREATE BIGWIGS WITH DEEPTOOLS
 	if (!params.egs) {
 		process create_fwd_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_fwd
 		val egs_size from egs_size_deeptools_fwd

 		output:
 		file("${id}.RPGCnorm.fwd.bigWig") into fwd_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --filterRNAstrand forward
 		"""
 	}

 	process create_rev_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_rev
 		val egs_size from egs_size_deeptools_rev

 		output:
 		file("${id}.RPGCnorm.rev.bigWig") into rev_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --filterRNAstrand reverse
 		"""
 	}
 }

 if (params.egs == '*') {
 		process create_fwd_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_fwd

 		output:
 		file("${id}.RPGCnorm.fwd.bigWig") into fwd_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${params.egs} --filterRNAstrand forward
 		"""
 	}

 	process create_rev_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_rev

 		output:
 		file("${id}.RPGCnorm.rev.bigWig") into rev_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${params.egs} --filterRNAstrand reverse
 		"""
 	}
 }

 	// STEP 6 QUALITY CONTROL WITH QORTS
 	process qorts {

 		publishDir "${params.outdir}/qc", mode: 'copy'

 		input:
 		file gtf_file
 		file fasta_file
 		set mergeid, id, file(bam), file(bam_index) from bam_qorts
 		set mergeid, id, file(read1) from qorts_trimmed_fastqs

 		output:
 		file 'QC/*' into qc_files

 		script:
 		"""
 		java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --randomSeed 111 --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id}  --singleEnded --rawfastq ${read1} ${bam} ${gtf_file} QC
 		"""
 	}

 	// STEP 7 PRESEQ
 	process preseq {

 		publishDir "${params.outdir}/preseq", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_preseq

 		output:
 		file("${id}.c_curve.txt") into preseq_results

 		script:
 		"""
 		preseq c_curve -l 9000000000 -B -o ${id}.c_curve.txt ${bam}
 		"""
 	}

 	// STEP 8 STRINGTIE
 	process stringtie {

 		publishDir "${params.outdir}/stringtie", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_stringtie
 		file gtf_file

 		output:
 		file "${id}.gene_abundance.txt"
 		file "${id}.cov_refs.gtf"
 		file "stringtie_log" into stringtie_log

 		script:
 		"""
 		stringtie ${bam} -G ${gtf_file} -A ${id}.gene_abundance.txt -C ${id}.cov_refs.gtf -e -b ${id}_ballgown -p ${params.threads} > stringtie_log
 		"""
 	}

 	// STEP 9 READ COUNTING WITH FEATURECOUNTS
 	if (params.strandInfo == 'unstranded') {

 		process featurecounts_unstranded {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_unstranded
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam}
 			featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 0 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	if (params.strandInfo == 'frFirstStrand') {

 		process featurecounts_frFirstStrand {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_frfirst
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam}
 			featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 2 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	if (params.strandInfo == 'frSecondStrand') {

 		process featurecounts_frSecondStrand {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_frsecond
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam}
 			featureCounts -a ${gtf_file} -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 1 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	// STEP 10 DGE with RUVSeq and EdgeR annd DESeq2
 	process dge {

 		publishDir "${params.outdir}/dge", mode: 'copy'

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

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()
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

 	} // closing bracket SE RNA

 	// PE RNA
 	if (params.mode == 'rna' && params.lib == 'p') {

 		if (!params.strandInfo) {
 	exit 1, "Please specify strand information. Available: unstranded, frFirstStrand, frSecondStrand. If you are unsure, run the pipeline using --strandInfo unstranded and --subsample and then look in your qc folder for information on the strandedness of your dataset."
 }

 		if (!params.expInfo) {
 	exit 1, "Please specify a experiment config file. Check README for exact information. Typically includes condition information for differential gene expression analysis."
 }

 exp_file = file(params.expInfo)

 	// Parse config file
 	fastqs = Channel
 	.from(config_file.readLines())
 	.map { line ->
 			list = line.split()
 			mergeid = list[0]
 			id = list[1]
 			path1 = file(list[2])
 			path2 = file(list[3])
 			[ mergeid, id, path1, path2 ]
 		}

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2) from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz") into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_fwd
 			egs_size_deeptools_rev
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'hisat2') {
	// Generate HISAT2 Index
 	process create_mapping_index {

 		publishDir "${params.outdir}/hisat2_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("*") into hisat22_index

 		script:
 		"""
 		hisat2-build -p ${params.threads} ${fasta_file} genome
 		"""
 	}
}

	if (params.aligner == 'star') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz") into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs, qorts_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 intronlen=${params.intronlen} maxindel=${params.maxindel} ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 // STEP 4 MAPPING WITH HISAT2
if (params.aligner == 'hisat2') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from bbmap_trimmed_fastqs
 		file("*") from hisat2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond
 		file("${id}.hs2.metricsFile.txt")

 		script:
 		"""
 		hisat2 -5 ${params.hs_trim5} -3 ${params.hs_trim3} --mp ${params.hs_mp} --sp ${params.hs_sp} --np ${params.hs_np} --rdg ${params.hs_rdg} --rfg ${params.hs_rfg} --pen-cansplice ${params.hs_pen_cansplice} --pen-noncansplice ${params.hs_pen_noncansplice} --min-intronlen ${params.hs_min_intronlen} --max-intronlen ${params.hs_max_intronlen} -k ${params.hs_k} --max-seeds ${params.hs_max_seeds} --met-file ${id}.hs2.metricsFile.txt -p ${params.threads} -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 if (params.aligner == 'star') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2) from bbmap_trimmed_fastqs
 		file("indexFiles/*") from star_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_deeptools_fwd, bam_deeptools_rev, bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond

 		script:
 		"""
 		STAR --genomeDir indexFiles --runThreadN ${params.threads} --readFilesIn ${read1} ${read2} --readFilesCommand zcat --clip3pNbases ${params.star_clip3pNbases} --clip5pNbases ${params.star_clip5pNbases} --outFileNamePrefix ${id} --outFilterMultipmapScoreRange ${params.star_outFilterMultimapScoreRange} --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} --outFilterScoreMin ${params.star_outFilterScoreMin} --alignEndsType ${params.star_alignEndsType} --winAnchorMultimapNmax ${params.star_winAnchorMultimapNmax} --outSAMtype BAM SortedByCoordinate --quantMode ${params.star_quantMode} --twopassMode ${params.star_twopassMode}
 		sambamba sort --tmpdir $baseDir -N -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}Aligned.sortedByCoord.out.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// STEP 5 CREATE BIGWIGS WITH DEEPTOOLS
 	process create_fwd_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_fwd
 		val egs_size from egs_size_deeptools_fwd

 		output:
 		file("${id}.RPGCnorm.fwd.bigWig") into fwd_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.fwd.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --filterRNAstrand forward
 		"""
 	}

 	process create_rev_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_deeptools_rev
 		val egs_size from egs_size_deeptools_rev

 		output:
 		file("${id}.RPGCnorm.rev.bigWig") into rev_bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.rev.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --filterRNAstrand reverse
 		"""
 	}

 	// STEP 6 QUALITY CONTROL WITH QORTS
 	process qorts {

 		publishDir "${params.outdir}/qc", mode: 'copy'

 		input:
 		file gtf_file
 		file fasta_file
 		set mergeid, id, file(bam), file(bam_index) from bam_qorts
 		set mergeid, id, file(read1), file(read2) from qorts_trimmed_fastqs

 		output:
 		file 'QC/*' into qc_files

 		script:
 		"""
 		java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --keepMultiMapped --title ${id} --randomSeed 111 --stopAfterNReads 5000000 --nameSorted --outfilePrefix ${id} --rawfastq ${read1},${read2} ${bam} ${gtf_file} QC
 		"""
 	}

 	// STEP 7 PRESEQ
 	process preseq {

 		publishDir "${params.outdir}/preseq", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_preseq

 		output:
 		file("${id}.c_curve.txt") into preseq_results

 		script:
 		"""
 		preseq c_curve -l 9000000000 -B -o ${id}.c_curve.txt -P ${bam}
 		"""
 	}

 	// STEP 8 STRINGTIE
 	process stringtie {

 		publishDir "${params.outdir}/stringtie", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), file(bam_index) from bam_stringtie
 		file gtf_file

 		output:
 		file "${id}.gene_abundance.txt"
 		file "${id}.cov_refs.gtf"
 		file "stringtie_log" into stringtie_log

 		script:
 		"""
 		stringtie ${bam} -G ${gtf_file} -A ${id}.gene_abundance.txt -C ${id}.cov_refs.gtf -e -b ${id}_ballgown -p ${params.threads} > stringtie_log
 		"""
 	}

 	// STEP 9 READ COUNTING WITH FEATURECOUNTS
 	if (params.strandInfo == 'unstranded') {

 		process featurecounts_unstranded {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_unstranded
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 0 ${bam}
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 0 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	if (params.strandInfo == 'frFirstStrand') {

 		process featurecounts_frFirstStrand {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_frfirst
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 2 ${bam}
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 2 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	if (params.strandInfo == 'frSecondStrand') {

 		process featurecounts_frSecondStrand {

 			publishDir "${params.outdir}/featurecounts", mode: 'copy'

 			input:
 			set mergeid, id, file(bam), file(bam_index) from bam_fc_frsecond
 			file gtf_file

 			output:
 			file "${id}_gene.featureCounts.txt" into geneCounts
 			file "${id}_gene.featureCounts.txt.summary" into featureCounts_logs
 			file "${id}_biotype_counts.txt" into featureCounts_biotype

 			script:
 			"""
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_id -o ${id}_gene.featureCounts.txt -s 1 ${bam}
 			featureCounts -a ${gtf_file} -T ${params.threads} -p -g gene_biotype -o ${id}_biotype.featureCounts.txt -s 1 ${bam}
 			cut -f 1,7 ${id}_biotype.featureCounts.txt > ${id}_biotype_counts.txt
 			"""

 		}
 	}

 	// STEP 10 DGE with RUVSeq and EdgeR annd DESeq2
 	process dge {

 		publishDir "${params.outdir}/dge", mode: 'copy'

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

 	// STEP 11 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()
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

 	} // closing bracket PE RNA

 	// SE ATAC
if (params.mode == 'atac' && params.lib == 's') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} out=${id}.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} out=${id}_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc ${read1}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 1 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno, np_NI_motifs

 		script:
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAM -g ${egs_size} -q ${params.qvalue} --nomodel --extsize 73 --shift -37 --seed 111 --broad --keep-dup all
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket SE atac

 	// PE ATAC
 if (params.mode == 'atac' && params.lib == 'p') {

 	// Parse config file
 	fastqs = Channel
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

 	// Subsample
 	if (params.subsample == true) {
 		process subsampling {

 			input:
 			set mergeid, id, file(read1), file(read2), controlid, mark from fastqs

 			output:
 			set mergeid, id, file("${id}_1.subsampled.fastq.gz"), file("${id}_2.subsampled.fastq.gz"), controlid, mark into subsampled_fastqs

 			script:
 			"""
 			reformat.sh in=${read1} in2=${read2} out=${id}_1.subsampled.fastq.gz out2=${id}_2.subsampled.fastq.gz sample=100000
 			"""
 		}

 		subsampled_fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}	else {
 		fastqs.into {
 			fastqs_fastqc
 			fastqs_bbduk
 		}
 	}

 	// Fetch chromosome sizes from FASTA
 	process fetch_chrom_sizes {

 		input:
 		file fasta_file

 		output:
 		file "chromSizes.txt" into chrom_sizes_WI, chrom_sizes_NI, chrom_sizes_epic

 		script:
 		"""
 		samtools faidx ${fasta_file}
 		awk -v OFS='\t' '{print \$1, \$2}' ${fasta_file}.fai > chromSizes.txt
 		"""
 	}

 	// Calculate effective genome size
 	process calculate_egs {

 		input:
 		file fasta_file

 		output:
 		file "egs_size.txt" into egs_file

 		script:
 		"""
 		epic-effective -r ${params.readLen} -n ${params.threads} -t $baseDir ${fasta_file} > egs_file.txt
 		Rscript $baseDir/bin/process_epic_effective_output.R egs_file.txt
 		"""
 	}

 	egs_file.map{ file ->
 		file.text.trim() } .set {
 			egs_size
 		}

 		egs_size.into {
 			egs_size_deeptools_NI
 			egs_size_macs_NI
 		}

 	// Generate BBMap Index
 	if (params.aligner == 'bbmap') {
 	process create_mapping_index {

 		publishDir "${params.outdir}/bbmap_index", mode: 'copy'

 		input:
 		file fasta_file

 		output:
 		file("ref/*") into bbmap_index

 		script:
 		"""
 		bbmap.sh ref=${fasta_file} usemodulo
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
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

	if (params.aligner == 'bwa') {
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

 	// STEP 1 PRE TRIM FASTQC
 	process fastqc_preTrim {

 		publishDir "${params.outdir}/preTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_fastqc

 		output:
 		file("*.zip") into pre_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 2 TRIMMING WITH BBDUK
 	process trimming {

 		publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqs_bbduk

 		output:
 		set mergeid, id, file("${id}_1_postTrimmed.fastq.gz"), file("${id}_2_postTrimmed.fastq.gz"), controlid, mark into bbmap_trimmed_fastqs, fastqc_trimmed_fastqs

 		script:
 		"""
 		bbduk.sh in=${read1} in2=${read2} out=${id}_1_postTrimmed.fastq.gz out2=${id}_2_postTrimmed.fastq.gz ref=$baseDir/adapters/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 		"""
 	}

 	// STEP 3 POST TRIM FASTQC
 	process fastqc_postTrim {

 		publishDir "${params.outdir}/postTrim_fastqc", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from fastqc_trimmed_fastqs

 		output:
 		file("*.zip") into post_fastqc_results

 		script:
 		"""
 		fastqc -t 2 ${read1} ${read2}
 		"""
 	}

 	// STEP 4 MAPPING WITH BBMAP
 	if (params.aligner == 'bbmap') {
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} in2=${read2} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 maxindel=1 ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	if (params.aligner == 'bowtie2') {
	// STEP 4 MAPPING WITH Bowtie2
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bowtie2_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		if (params.local == true) {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} --local -x genome -1 ${read1} -2 ${read2} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		} else {
 		"""
 		bowtie2 -q -D ${params.bt2_D} -R ${params.bt2_R} -N ${params.bt2_N} -L ${params.bt2_L} -i ${params.bt2_i} -5 ${params.bt2_trim5} -3 ${params.bt2_trim3} -k ${params.bt2_k} -p ${params.threads} -x genome -U ${read1} -S ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 		}
 	}
 }

 	if (params.aligner == 'bwa') {
 	// STEP 4 MAPPING WITH BWA
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1), file(read2), controlid, mark from bbmap_trimmed_fastqs
 		file("*") from bwa_index

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), controlid, mark, file("${id}.sorted.mapped.bam.bai") into bam_grouping, bam_spp

 		script:
 		"""
 		bwa mem -t ${params.threads} -k ${params.bwa_k} -w ${params.bwa_w} -d ${params.bwa_d} -r ${params.bwa_r} -c ${params.bwa_c} -A ${params.bwa_A} -B ${params.bwa_B} -O ${params.bwa_O} -E ${params.bwa_E} -L ${params.bwa_L} -U ${params.bwa_U} -T ${params.bwa_T} -M genome ${read1} ${read2} > ${id}.mapped.sam
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.sam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}
 }

 	// SEPARATE CHIP AND INPUT FILES
 	//treat = Channel.create()
 	//control = Channel.create()
 	//bam_grouping.choice(treat, control) {
 	//	it[4] == 'input' ? 1 : 0
 	//}

 	// STEP 5 ESTIMATE FRAGMENT SIZE SPP
 	process estimate_fragment_size {

 		publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), controlid, mark, file(bam_index) from bam_grouping

 		output:
 		set mergeid, id, file("${id}.params.out") into modelParams
 		set mergeid, id, file(bam), controlid, file("${id}.params.out"), mark, file(bam_index) into modelBams

 		script:
 		"""
 		Rscript $baseDir/bin/run_spp.R -c=${bam} -rf -out=${id}.params.out -savp=${id}.pdf -p=${params.threads} -tmpdir=$baseDir
 		"""
 	}

 	(bams) = modelBams.map { mergeid, id, bam, controlid, paramFile, mark, bam_index ->
 		fragLen = paramFile.text.split()[2].split(',')[0]
 		[ mergeid, id, bam, controlid, mark, fragLen, bam_index ]
 	}.into(1)

 	bams.into {
 		bams_bigwigs
 		bams_macs
 	}

 	// STEP 6 GENERATE RPGC NORMALIZED COVERAGE TRACKS WITH NO INPUT
 	process create_coverage_tracks {

 		publishDir "${params.outdir}/tracks", mode: 'copy'

 		input:
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_bigwigs
 		val egs_size from egs_size_deeptools_NI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 1 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReads
 		"""
 	}

 	// STEP 7 NARROW PEAK CALLING WITH MACS2 NO INPUT
 	process call_binding_sites {

 		publishDir "${params.outdir}/peaks", mode: 'copy'

 		input:
 		file chromSizes from chrom_sizes_NI.val
 		set mergeid, id, file(bam), mark, fragLen, file(bam_index) from bams_macs
 		val egs_size from egs_size_macs_NI

 		output:
 		set mergeid, id, file("${id}/${id}_peaks.narrowPeak"), mark, fragLen into np_NI_anno, np_NI_motifs

 		script:
 		"""
 		macs2 callpeak -t ${bam} -n ${id} --outdir ${id} -f BAMPE -g ${egs_size} --tempdir $baseDir -q ${params.qvalue} --nomodel --extsize 73 --shift 37 --broad --seed 111
 		"""
 	}

 	// STEP 9 ANNOTATE PEAKS USING HOMER
 	process annotate_binding_sites {

 		publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

 		input:
 		file fasta_file
 		file gtf_file
 		set mergeid, id, file(np), mark, fragLen from np_NI_anno

 		output:
 		file("${id}_narrowPeaks.annotated.txt")

 		script:
 		"""
 		perl $baseDir/bin/homer/bin/annotatePeaks.pl ${np} ${fasta_file} -gtf ${gtf_file} > ${id}_narrowPeaks.annotated.txt
 		"""
 	}

 	// STEP 10 MULTIQC
 	process multiqc {

 		publishDir "${params.outdir}/multiqc", mode: 'copy'

 		input:
 		file ('fastqc/*') from post_fastqc_results.flatten().toList()
 		file ('fastqc/*') from pre_fastqc_results.flatten().toList()

 		output:
 		file "*multiqc_report.html"
 		file "*multiqc_data"

 		script:
 		"""
 		multiqc -f .
 		"""
 	}

 	} // closing bracket PE atac

 	// ANALYSIS MODE
 	if (params.mode == 'analysis' && params.analysis == 'predictEnhancers') {

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

			publishDir "${params.outdir}/predicted_enhancers", mode: 'copy'

			input:
			file(enhancer_file) from enhancers

			output:
			file("predicted_enhancers.txt")

			script:
			"""
			sort-bed ${enhancer_file} | bedops --merge - | awk -v OFS='\t' '{print \$1, \$2, \$3, "Enhancer.Prediction." NR}' > predicted_enhancers.txt
			"""
		}

 	} // closing bracket ANALYSIS mode

 // ON COMPLETION
 workflow.onComplete {
 	println ""
 	println "Workflow completed on: $workflow.complete"
 	println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
 	println "Workflow Duration: $workflow.duration"
 	println ""
 	println "Submit issues/requests on GitHub or e-mail cag104@ucsd.edu"
 	println ""
 }
