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
 //params.binSize = 10
 //params.smoothLen = 50
 //params.aligner = 'bbmap'

 // PRINT HELP
 if (params.help) {
 	log.info ''
 	log.info '~ C I P H E R ~ Version 1.0.0'
 	log.info '**************************************************************'
 	log.info ''
 	log.info 'REQUIRED PARAMETERS:'
 	log.info '--mode				Choose from available: chip, rna, gro, mnase, dnase, analysis.'
 	log.info '--config			Configuration file with sample information. Check README for more information.'
 	log.info '--fasta				Reference genome in FASTA format.'
 	log.info '--gtf				Reference genome in GTF format.'
 	log.info '--lib				Library information. s for single-stranded data, p for pair-ended data.'
 	log.info '--readLen			The length of your reads.'
 	log.info ''
 	log.info 'RNA-seq ONLY:'
 	log.info '--strandInfo		Strandedness information. Choose from "unstranded", "frFirstStrand", "frSecondStrand".'
 	log.info ''
 	log.info 'OPTIONAL PARAMETERS:'
 	log.info '--threads			Number of threads. (Default: 1)'
 	log.info '--minid				Minimum alignment identity to look for. Higher is faster and less sensitive. (Default: 0.76)'
 	log.info '--qvalue			Minimum FDR cutoff for peak detection. (Default: 0.05)'
 	log.info '--epic_w			Size of the windows used to scan the genome for peak detection in EPIC. (Default: 200)'
 	log.info '--epic_g			A multiple of epic_w used to determine the gap size in EPIC. (Default: 3)'
 	log.info '--maxindel		Maximum indel length searched. 200k recommended for vertebrate genomes. (Default: 200k)'
 	log.info '--intronlen		Maximum intron length. 20 recommended for vertebrate genomes. (Default: 20)'
 	log.info '--outdir			Name of output directory. (Default: results)'
 	log.info ''
 	log.info '--subsample			Set this flag to subsample reads for testing.'
 	log.info ''
 	log.info '**************************************************************'
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
 log.info ''
 log.info "subsample:	${params.subsample}"
 log.info ''
 log.info '**********************************************'
 log.info ''

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
 			egs_size_deeptools_WI
 			egs_size_deeptools_NI
 			egs_size_macs_WI
 			egs_size_macs_NI
 			egs_size_epic_WI
 		}

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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_bigwigs
 		val egs_size from egs_size_deeptools_WI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_WI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReadss
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_macs
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_epic
 		val egs_size from egs_size_epic_WI

 		output:
 		set mergeid, id, file("${id}_epic.broadPeak"), mark, fragLen into bp_WI_anno

 		script:
 		"""
 		bedtools bamtobed -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_size} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.qvalue} -cs ${chromSizes} > ${id}_epic.broadPeak
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_chipqc
 		set mergeid, id, file(np), mark, fragLen np_WI_chipqc

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
 			egs_size_epic_WI
 		}

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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_bigwigs
 		val egs_size from egs_size_deeptools_WI

 		output:
 		file("${id}.RPGCnorm.bigWig") into bigwigs_WI

 		script:
 		"""
 		bamCoverage -b ${bam} -o ${id}.RPGCnorm.bigWig -of bigwig -bs 10 -p ${params.threads} --normalizeTo1x ${egs_size} --smoothLength 50 -e ${fragLen} --ignoreDuplicates --centerReadss
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_macs
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_epic
 		val egs_size from egs_size_epic_WI

 		output:
 		set mergeid, id, file("${id}_epic.broadPeak"), mark, fragLen into bp_WI_anno

 		script:
 		"""
 		bedtools bamtobed -bedpe -i ${bam} > ${id}_treatment.bed
 		bedtools bamtobed -bedpe -i ${control} > ${id}_control.bed
 		epic --treatment ${id}_treatment.bed --control ${id}_control.bed -cpu ${params.threads} -egs ${egs_size} --fragment-size ${fragLen} -w ${params.epic_w} -g ${params.epic_g} -fdr ${params.qvalue} -cs ${chromSizes} --pair-end > ${id}_epic.bed
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
 		set mergeid, id, file(bam), file(control), mark, fragLen, file(bam_index), file(control_index) from bamsWI_chipqc
 		set mergeid, id, file(np), mark, fragLen np_WI_chipqc

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
 		python $baseDir/bin/danpos/danpos.py dpos ${bed} -m 0 -o ${id} --frsz ${fragLen} -jw 40 -jd 150
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
 			egs_size_deeptools_fwd
 			egs_size_deeptools_rev
 		}

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
 		python $baseDir/bin/danpos/danpos.py dpos ${bed} -m 1 -o ${id} --frsz ${fragLen} -jw 40 -jd 150
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
 	exit 1, "Please specify strand information. Available: unstranded, frFirstStrand, frSecondStrand. If you are unsure, run the pipeline using --strandInfo unstranded and --subsample and then look in your qc folder for information."
 }

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
 	process mapping {

 		publishDir "${params.outdir}/alignments", mode: 'copy'

 		input:
 		set mergeid, id, file(read1) from bbmap_trimmed_fastqs
 		file("ref/*") from bbmap_index.first()

 		output:
 		set mergeid, id, file("${id}.sorted.mapped.bam"), file("${id}.sorted.mapped.bam.bai") into bam_qorts, bam_preseq, bam_stringtie, bam_fc_unstranded, bam_fc_frfirst, bam_fc_frsecond
 		file("${id}.alignmentReport.txt")
 		file("${id}.unmapped.bam") into unmapped_bams

 		script:
 		"""
 		bbmap.sh in=${read1} outm=${id}.mapped.bam outu=${id}.unmapped.bam keepnames=t trd sam=1.3 intronlen=${params.intronlen} maxindel=${params.maxindel} ambig=random statsfile=${id}.alignmentReport.txt minid=${params.minid} usemodulo
 		sambamba sort --tmpdir $baseDir -t ${params.threads} -o ${id}.sorted.mapped.bam ${id}.mapped.bam
 		sambamba index -t ${params.threads} ${id}.sorted.mapped.bam
 		"""
 	}

 	// STEP 5 QUALITY CONTROL WITH QORTS
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
 		java -jar $baseDir/bin/QoRTs.jar QC --generatePlots --genomeFA ${fasta_file} --unstranded --title ${id} --randomseed 111 --outfilePrefix ${id}  --singleEnded --rawfastq ${read1} ${bam} ${gtf_file} QC
 		"""
 	}

 	// STEP 6 PRESEQ
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

 	// STEP 7 STRINGTIE
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

 	// STEP 8 READ COUNTING WITH FEATURECOUNTS
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

 		process featurecounts_unstranded {

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

 		process featurecounts_unstranded {

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

 	// STEP 9 DGE with RUVSeq and EdgeR annd DESeq2

 	} // closing bracket SE RNA

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
