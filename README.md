# bare_bcbio
a container that runs some bcbio pipelines with no wrapper script and minimal input

This container implements a somatic variant calling and germline variant calling pipeline run internally by bcbio. See the somatic variant calling pipeline described in the bcbio documentation at:
http://bcbio-nextgen.readthedocs.io/en/latest/contents/testing.html#cancer-tumor-normal

IF CWLTOOL OR DOCKSTORE IS USED TO RUN THE CONTAINER SET THE HOST ENVIRONMENT VARIABLE 'TMPDIR' TO A DIRECTORY WITH ENOUGH SPACE TO ACCOMODATE VERY LARGE FILES


Sample command line to run somatic variant calling using the docker run command: 
```
    docker run -it 
    -v /data/input/:/data/input/ 
    -v /data/ops/:/data/ops/ 
    -v /home/ubuntu/GATK/GenomeAnalysisTK-3.6.tar.bz2:/home/ubuntu/GATK/GenomeAnalysisTK-3.6.tar.bz2 
    -v /reference/data/:/reference/data/ 
    -v $(pwd):$(pwd) 
    -w $(pwd) 
    -u 1001:1001 
    quay.io/wshands/bare_bcbio:latest UCSCbcbioTool.py 
    -n /data/input/SPCG-HB453_13G_normal_R1.fastq.gz 
    -n /data/input/SPCG-HB453_13G_normal_R2.fastq.gz 
    -t /data/input/SPCG-HB453_8R_tumor_R1.fastq.gz 
    -t /data/input/SPCG-HB453_8R_tumor_R2.fastq.gz 
    -g /home/ubuntu/GATK/GenomeAnalysisTK-3.6.tar.bz2 
    -b /data/ops/NGv3.bed 
    -W somatic-variant-calling 
    -o $(pwd)/final 
    -c 24 
    -d /reference/data/ 
    > out
```
The work directory, where temp files will be created, specified by -w, must be "mirror mounted". This means that a -v switch must be used to provide the container access to the working directory and the directory must be the same on both sides of the colon in Docker's -v command.

Additionally, the path to all input files or data locations must be "mirror mounted" with a -v switch.

Use the Docker '-u' switch with your user id and group id to enable the container to run as the current user. One way to do this is to use the argument $(grep $USER /etc/passwd | cut -d: -f3,4) like this:
```
    -u $(grep $USER /etc/passwd | cut -d: -f3,4)
```
If the current user also is the owner of the current working directory this will allow the container to write to the working directory (which is set using the '-w' switch).

You may do structural variant calling at the same time as either somatic or germ-line small variant calling by specifying both on the command line with the -W switch, e.g.:
```
    ...
    -W somatic-variant-calling 
    -W structural-variant-calling 
    ...
```
You may point to a directory where bcbio reference data is located, or an empty directory where it will be downloaded, with the -d switch, or you may point to a tar ball of the bcbio reference data which will be un-tarred in the current working directory with the -t switch, e.g.
```
    -d /reference/data/empty/dir/
```
or
```
    -d /reference/data/already/downloaded/
```
or
```
    -t /reference/data/tar/ball/file.tar
```
Help output for python script inside container:
```
docker run -it quay.io/wshands/bare_bcbio UCSCbcbioTool.py --help
usage: UCSCbcbioTool.py [-h] [-t TUMOR_INPUT] [-n NORMAL_GERMLINE_INPUT]
                        [-c NUM_CORES] -g GATK_FILE
                        [-W {somatic-variant-calling,germline-variant-calling}]
                        [-b WES_BED_FILE] [-s] [-r RUN_NAME]
                        [-d DATA_DIR | -f DATA_FILE] [-o OUTPUT_DIR]

Script to run the bcbio Docker container

optional arguments:
  -h, --help            show this help message and exit
  -t TUMOR_INPUT, --tumor_input TUMOR_INPUT
                        Input file(s) for processing. Multiple path/file names
                        can be provided separated by spaces. These could
                        correspond fastq files for tumor paired end reads.
  -n NORMAL_GERMLINE_INPUT, --normal_germline_input NORMAL_GERMLINE_INPUT
                        Input file(s) for processing. Multiple path/file names
                        can be provided separated by spaces. These could
                        correspond fastq files for normal paired end reads.
  -c NUM_CORES, --num_cores NUM_CORES
                        number of cores to use for processing.
  -g GATK_FILE, --GATK_file GATK_FILE
                        Path to GATK file, e.g.
                        /path/to/GenomeAnalysisTK.tar.bz2.
  -W {somatic-variant-calling,germline-variant-calling}, --workflow {somatic-variant-calling,germline-variant-calling}
                        The name of the workflow to run. Default is somatic
                        variant calling.
  -b WES_BED_FILE, --WES_bed_file WES_BED_FILE
                        Input BED file for WES variant calling. Include if
                        whole exome (WES) calling is desired. Default is whole
                        genome (WGS) calling
  -s, --include_sv      Include structural variant calling in workflow.
                        Default is False
  -r RUN_NAME, --run_name RUN_NAME
                        Name to use for batching samples and intermediate file
                        names.
  -d DATA_DIR, --data_dir DATA_DIR
                        Directory where reference genome files are or shouldbe
                        downloaded. Default is to download the reference files
                        to the current working directory.
  -f DATA_FILE, --data_file DATA_FILE
                        Path to reference genomes tar file. Default is to
                        download the reference files to the current working
                        directory.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory where output files should be written.
                        Default is <cwd>/final
```
