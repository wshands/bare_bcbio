#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "UCSCbcbioTool"
label: "unwrapped bcbio nextgen tool"
cwlVersion: v1.0 
doc: |
    ![build_status](https://quay.io/wshands/bare_bcbio/status)
    A Docker container for the bcbio command. 
    See the bcbio (http://bcbio-nextgen.readthedocs.io/en/latest/index.html) 
    website for more information.
    ```
     Usage:
    # fetch CWL
    $> dockstore tool cwl --entry quay.io/wshands/bare_bcbio/bare_bcbio_germline_variant_calling:1.0.3 > germline_variant_calling.cwl
    # make a runtime JSON template and edit it
    $> dockstore tool convert cwl2json --cwl germline_variant_calling.cwl > germline_variant_calling.json
    # run it locally with the Dockstore CLI
    $> dockstore tool launch --debug  --entry quay.io/wshands/bare_bcbio/bare_bcbio_germline_variant_calling:1.0.3  --json germline_variant_calling.json
    ```

dct:creator:
  foaf:name: Walt Shands
  foaf:mbox: "jshands@ucsc.edu"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wshands/bare_bcbio:1.0.3"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000
    description: "the process requires at least 4G of RAM"

inputs:
  normal_germline_input:
    type:
      type: array
      items: File
      doc: "Somatic calling normal input file(s) for processing"
      inputBinding:
        prefix: --normal_germline_input

  num_cores:
    type: int?
    default: 16
    doc: "The number of cores bcbio should use for processing"
    inputBinding:
      prefix: --num_cores

  GATK_file:
    type: File
    doc: "Path to GATK file, e.g. /path/to/GenomeAnalysisTK.tar.bz2"
    inputBinding:
      prefix: --GATK_file

  workflow:
    type: string?
    doc: "Workflow to run; [somatic_variant_calling, germline_variant_calling]"
    default: "germline-variant-calling"
    inputBinding:
      prefix: --workflow

  WES_bed_file:
    type: File?
    doc: "Input BED file. Include if whole exome (WES) calling is desired."
    inputBinding:
      prefix: --WES_bed_file

  include_sv:
    type: boolean?
    default: false
    doc: "Include structural variant calling"
    inputBinding:
      prefix: --include_sv

  run_name:
    type: string?
    default: "Current_run"
    doc: "Name to use for batching samples and intermediate file names."
    inputBinding:
      prefix: --run_name
 
  data_file:
    type: File?
    doc: "Input tar file of bcbio reference genomes directory"
    inputBinding:
      prefix: --data_file

outputs:
  output_files:
    type:
      type: array
      items: File
    outputBinding:
      # should be put in the working directory
       glob: ./*
    doc: "Result files from the work flow"


baseCommand: ["UCSCbcbioTool.py"]
