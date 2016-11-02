#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "UCSCbcbioTool"
label: "bcbio nextgen tool"
cwlVersion: v1.0 
doc: |
    ![build_status](https://quay.io/repository/collaboratory/dockstore-tool-bamstats/status)
    A Docker container for the bcbio command. 
    See the bcbio (http://bcbio-nextgen.readthedocs.io/en/latest/index.html) 
    website for more information.
    ```
    Usage:
    # fetch CWL
    $> dockstore tool cwl --entry quay.io/collaboratory/UCSCbcbioTool > UCSCbcbioTool.cwl
    # make a runtime JSON template and edit it
    $> dockstore tool convert cwl2json --cwl UCSCbcbioTool.cwl > UCSCbcbioTool.json
    # run it locally with the Dockstore CLI
    $> dockstore tool launch --entry quay.io/collaboratory/UCSCbcbioTool:latest \
        --json UCSCbcbioTool.json
    ```

#dct:creator:
#  "@id": "jshands@ucsc.edu"
#  foaf:name: Walt Shands
#  foaf:mbox: "jshands@ucsc.edu"

requirements:
  - class: DockerRequirement
    dockerPull: "ucscbcbiotool:latest"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000
    description: "the process requires at least 4G of RAM"

inputs:
  sample_files:
    type:
      type: array
      items: File
      doc: "Germline calling input file(s) for processing"
      inputBinding:
        prefix: --sample_files

  num_cores:
    type: int?
    default: 16
    doc: "The number of cores bcbio should use for processing"
    inputBinding:
      prefix: --num_cores

  GATK_file:
    type: File?
    doc: "Path to GATK file, e.g. /path/to/GenomeAnalysisTK.tar.bz2"
    inputBinding:
      prefix: --GATK_file

  workflow:
    type: string
    doc: "Workflow to run; [germline_variant_calling, structural_variant_calling]"
    inputBinding:
      prefix: --workflow

  bed_file:
    type: File
    doc: "Input BED file"
    inputBinding:
      prefix: --bed_file

  data_file:
    type: File?
    doc: "Input tar file of bcbio reference genomes directory"
    inputBinding:
      prefix: --data_file

utputs:
  output_files:
    type:
      type: array
      items: File
    outputBinding:
      # should be put in the working directory
       glob: ./*
    doc: "Result files from the work flow"


baseCommand: ["UCSCbcbioTool.py"]
