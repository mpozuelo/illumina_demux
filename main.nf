#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/illumina_demux
========================================================================================
 mpozuelo/illumina_demux Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/mpozuelo/illumina_demux
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/illumina_demux --input '*.txt' -profile docker

    Mandatory arguments:
      --run [file]                  Illumina run we are demultiplexing
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.
      --cluster_path                Path in the cluster for the genomics unit storage (Default: "/datos/ngs/dato-activo/data/")

    QC:
      --skipQC                      Skip all QC steps apart from MultiQC
      --skipFastQC                  Skip FastQC

    Trimming:
      --complete                    Makes the STAR with the whole sample (for example for culture cell CIMA contamination)

    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }

run = params.run
sequencer = params.sequencer
cluster_path = params.cluster_path


runDir = file("${cluster_path}/data/01_bcl/Illumina/$sequencer/$run", checkIfExists: true)
ch_samplesheet = file("${runDir}/SampleSheet.csv", checkIfExists: true)
protocol = params.protocol

// Validate inputs
cluster_path = params.cluster_path


if (!params.outdir) {
  params.outdir = params.run
}




// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Name'] = custom_runName ?: workflow.runName
summary['Run'] = params.run
summary['Sequencer'] = params.sequencer
summary['Technology'] = "Illumina"
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-illumina_demux-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/illumina_demux Workflow Summary'
    section_href: 'https://github.com/mpozuelo/illumina_demux'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
Get the ONLY the samplesheet for getting information per sample for the Trimming
Also get information about the cycles for the demultiplexing
*/

process parse_samplesheet {
  tag "$samplesheet"
  label 'process_low'
  publishDir "${cluster_path}/data/03_intermediate/Illumina/${sequencer}/${run}/SampleSheet/", mode: 'copy',
  saveAs: { filename ->
    filename.startsWith("samplesheet_demux") ? filename : null
  }

  input:
  path samplesheet from ch_samplesheet

  output:
  path "samplesheet_demux*" into ch_samples_info
  //path "SampleSheet.csv" into ch_samplesheet_bcl2fastq2
  path ("*.txt") into ch_demux_parameters

  script:
  """
  echo \$(grep -A 1 'Reads' SampleSheet.csv | tail -n 1 | tr -dc '0-9') > cycles1.txt
  echo \$(grep -A 2 'Reads' SampleSheet.csv | tail -n 1 | tr -dc '0-9') > cycles2.txt
  echo \$(grep -A 3 'Reads' SampleSheet.csv | tail -n 1 | tr -dc '0-9') > cycles3.txt
  echo \$(grep -A 4 'Reads' SampleSheet.csv | tail -n 1 | tr -dc '0-9') > cycles4.txt
  grep -A 99999999999 'Data' SampleSheet.csv | grep -v Data > samplesheet_demux_run${run}.csv
  """
}


/*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing
 It contains the following columns:
0- Sample_ID (mandatory for Illumina)
1- Sample_Name (mandatory for Illumina)
2- index (mandatory for Illumina)
3- index2 (mandatory for Illumina)
4- SampleProject
5- SampleSource
6- Organism
7- Method
8- Library
9- Platform
10- Run
11- Date
12- User
13- Coverage
*/

/*
Channel
  .from( ch_samples_info )
  .splitCsv(header:false, sep:',')
  .map { it = ["${it[1]}", "${it[9]}", "${it[11]}",
  ]}
  .set { ch_fastqc }
  */



/*
 * STEP 1 - Demultiplex using bcl2fastq2
 */

/*
***For Illumina we will demultiplex with bcl2fastq2 software
*/

process demux {
  container 'mpozuelo/illuminademux:bcl2fastq'
  tag "$run"
  label 'process_high'
  publishDir "${cluster_path}/data/04_pfastq/Illumina/${sequencer}/${run}/", mode: 'copy'
  /*saveAs: { filename ->
    filename.endsWith(".fastq.gz") ? filename : "logs/$filename"
  }
  */

  input:
  path samplesheet from ch_samplesheet
  path(cycles) from ch_demux_parameters

  output:
  file "*.fastq.gz" into ch_fastqc
  file "Reports"
  file "Stats"
  file "*.log"

  script:
  info = "${run}.dmux.log 2>&1"

  """
  cycles1=\$(cat ${cycles[0]})
  cycles2=\$(cat ${cycles[1]})
  cycles3=\$(cat ${cycles[2]})
  cycles4=\$(cat ${cycles[3]})

  if [ protocol == "" ]
  then
  bases_mask=\$(printf "Y%s,I%s,I%s,Y%s" "\$cycles1" "\$cycles2" "\$cycles3" "\$cycles4")
  let minlength=\$cycles1-\$cycles2
  let short_adapter_read=\$cycles2-1
  elif [ protocol == "mcSCRBseq" ] 
  then
  let read1=\$cycles1-\$cycles2
  bases_mask=\$(printf "I%sY%s,I%s,N%s,Y%s" "\$cycles2" "\$read1" "\$cycles2" "\$cycles3" "\$cycles4")
  elif [ protocol == "marseq" ]
  then
  let i2=\$cycles3-\$cycles2
  let read2=\$cycles3-\$i2
  bases_mask=\$(printf "Y%s,I%s,I%sY%s" "\$cycles1" "\$cycles2" "\$i2" "\$read2")
  fi



  bcl2fastq \\
    --runfolder-dir ${runDir} \\
    --output-dir  ./  \\
    --use-bases-mask \$bases_mask \\
    --create-fastq-for-index-reads \\
    --sample-sheet $samplesheet \\
    --minimum-trimmed-read-length \$minlength \\
    --mask-short-adapter-read \$short_adapter_read \\
    --no-lane-splitting \\
    --barcode-mismatches 1 \\
    -r 8 \\
    -p 10 \\
    -w 10 \\
    -l INFO >> $info
  """
}
/*
fqname_fqfile_ch = ch_fastqc.map { fqFile -> [fqFile.getParent().getName(), fqFile ] }
Channel
  .from( ch_samples_info )
  .splitCsv(header:false, sep:',')
  .map { fqFile -> ["${fqFile[1]}", "${fqFile[4]}" ] }
  .set { ch_project }
//h_project = ch_samplesheet.map { fqFile -> ["${fqFile[1]}", "${fqFile[4]}" ] }
ch_fastqc_all = Channel.empty()
ch_fastqc_all = ch_fastqc_all.mix(fqname_fqfile_ch, ch_project)


/*
 * FastQC
 */
/*
 if (!params.skipQC || !params.skipFastQC) {

   process fastqc {
     tag "$sample"
     label 'process_medium'
     publishDir "${cluster_path}/data/04_rfastq/Illumina/${sequencer}/${run}/", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".zip") ? "zips/$filename" : "html/$filename"
     }

     input:
     set path(reads) from ch_fastqc

     output:
     path("*_fastqc.{zip,html}")

     script:
     """
     fastqc --quiet --threads $task.cpus $reads
     """
   }
 } else {
   fastqc_results = Channel.empty()
 }


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/illumina_demux] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/illumina_demux] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/illumina_demux]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/illumina_demux]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/illumina_demux v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
