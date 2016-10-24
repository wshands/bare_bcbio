#!/usr/bin/env python 
from __future__ import print_function, division 
""" 
    author Walt Shands 
    jshands@ucsc.com 

This script runs a cancer variant and or structural variant calling pipeline 
using the publicly available bcbio/bcbio Docker container.
See this document for more information on bcbio: 
        http://bcbio-nextgen.readthedocs.io/en/latest/index.html

It will download reference data if the data directory provided to it is empty. 
It will installthe GATK tools if the GATK switch and arguent is specified on the
 command line.

Inputs:
    The tumor fastq or BAM files: These should be separated by a space. Required.
    The normalfastq or BAM files: These should be separated by a space. Required.
    Total available cores. This tells bcbio how many total cores to use. 
               http://bcbio-nextgen.readthedocs.io/en/latest/contents/parallel.html
    The path and file name of the GATK tools. If provided the GATK tools will be installed
          in the container for use in the workflow. (required????) 
    The name of the workflow to run. If not provided cancer variant calling only is run.
    The path to and name of the BED file. Required.
    The path of the output directory. Required.     

"""

import sys, argparse, os
import collections
import time
import string
import subprocess
import tempfile
import multiprocessing
import getpass, grp, pwd            #for setting up correct user and group	

 
def parse_arguments():
    """
    Parse Command Line
    """
    parser = argparse.ArgumentParser(description='Script to run the bcbio Docker container')

#    parser.add_argument('-s', '--sample_files', type=str, action='append',
#                         help="Input file(s) for processing. Multiple path/file"  
#                         "names can be provided separated by spaces. These could" 
#                         "correspond fastq files for paired end reads." )

    parser.add_argument('-t', '--tumor_sample_files', type=str, required=True,
                         action='append',
                         help='Input file(s) for processing. Multiple path/file \
                         names can be provided separated by spaces. These could \
                         correspond fastq files for tumor paired end reads.' )

    parser.add_argument('-n', '--normal_sample_files', type=str, required=True,
                         action='append', help='Input file(s) for processing. \
                         Multiple path/file names can be provided separated by \
                         spaces. These could correspond fastq files for normal \
                         paired end reads.' )

    parser.add_argument('-c', '--num_cores',  type=int, default=16,  
                         help="number of cores to use for processing.") 

    parser.add_argument('-g', '--GATK_file', type=str, help='Path to GATK file,' \
                        'e.g. /path/to/GenomeAnalysisTK.tar.bz2.')
 
    parser.add_argument('-W','--workflow', choices=['cancer_variant_calling',
 #                         'germline-variant-calling','structural_variant_calling'], 
                         'structural_variant_calling'], 
                         default=['cancer_variant_calling'], 
                         action='append', 
                         help="The name of the workflow to run. Default is cancer" \
                         " variant calling.")

    parser.add_argument('-b', '--bed_file', type=str, required=True, 
                          help='Input BED file.')

    parser.add_argument('-o', '--output_dir', type=str, required=True, 
                       help='Directory where output files should be written.' )

    options = parser.parse_args()

    return (options)


def get_bcbio_system_template():
    """
    This templated describes the system resources to the bcbio Docker container
    so that it can allocate the correct resources to each thread.
    """
    bcbio_system_template = string.Template("""
#
resources:
  default:
    cores: $core_count
    jvm_opts:
    - -Xms750m
    - -Xmx3500m
    memory: 3G
  dexseq:
    memory: 10g
  express:
    memory: 8g
  gatk:
    jvm_opts:
    - -Xms500m
    - -Xmx3500m
  macs2:
    memory: 8g
  qualimap:
    memory: 4g
  seqcluster:
    memory: 8g
  snpeff:
    jvm_opts:
    - -Xms750m
    - -Xmx4g
galaxy_config: /mnt/biodata/galaxy/
    """)

    return bcbio_system_template



def get_cancer_variant_template():
    """
    This is the template that describes to bcbio the tools to used in the cancer
    variant calling workflow and if structural variant calling is also requested
    which tools to use. The input tumor and normal path and file names are 
    inserted into the template as is the path and file name of the BED file,
    the structural variant calling tools to use if structural variant calling is
    requested and the output directory.
    """
    cancer_variant_template = string.Template("""
# Cancer tumor/normal calling evaluation using synthetic dataset 3
# from the ICGC-TCGA DREAM challenge:
# https://www.synapse.org/#!Synapse:syn312572/wiki/62018
---
details:
- algorithm:
    aligner: bwa
    align_split_size: 5000000
    nomap_split_targets: 100
    mark_duplicates: true
    recalibrate: true
    realign: true
    remove_lcr: true
    platform: illumina
    quality_format: standard
    variantcaller: [mutect2, freebayes, vardict, varscan]
    indelcaller: false
    ensemble:
      numpass: 2
    variant_regions: $bed_file
#    variant_regions: /mnt/cancer-dream-syn3/input/NGv3.bed
    # svcaller: [cnvkit, lumpy, delly]
    $svcaller_info
    # coverage_interval: amplicon
  analysis: variant2
  description: syn3-normal

  # The YAML below should look like the following after substitution:
  # files: [ <path and file name of normal paired end reads>, <path and file name of other end of normal paired end reads> ]
  files: [$normal_sample_files]

  #files: /mnt/cancer-dream-syn3/input/synthetic.challenge.set3.normal.bam
#  files:
#    - /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_normal_NGv3_1.fq.gz
#    - /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_normal_NGv3_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: syn3
    phenotype: normal
- algorithm:
    aligner: bwa
    align_split_size: 5000000
    nomap_split_targets: 100
    mark_duplicates: true
    recalibrate: true
    realign: true
    remove_lcr: true
    platform: illumina
    quality_format: standard
    variantcaller: [mutect2, freebayes, vardict, varscan]
    indelcaller: false
    ensemble:
      numpass: 2
    variant_regions: $bed_file
#    variant_regions: /mnt/cancer-dream-syn3/input/NGv3.bed
#    validate: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf.gz
#    validate_regions: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth_regions.bed
    # svcaller: [cnvkit, lumpy, delly]
    $svcaller_info
    # coverage_interval: amplicon
  #   svvalidate:
  #     DEL: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_DEL.bed
  #     DUP: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_DUP.bed
  #     INS: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_INS.bed
  #     INV: /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_INV.bed
  analysis: variant2
  description: syn3-tumor


  # The YAML below should look like the following after substitution:
  # files: [ <path and file name of tumor paired end reads>, <path and file name of other end of tumor paired end reads> ]
  files: [$tumor_sample_files]

  #files: /mnt/cancer-dream-syn3/input/synthetic.challenge.set3.tumor.bam
#  files:
#    - /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_NGv3_1.fq.gz
#    - /mnt/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_NGv3_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: syn3
    phenotype: tumor
fc_date: '2014-08-13'
fc_name: dream-syn3
upload:
#  dir: /mnt/cancer-dream-syn3/final
  dir: $output_dir
    """)

    return cancer_variant_template

def get_user_group_for_a_path(path_to_check):
    """
    Get the effective user and group names and ids of the person
    who owns a particular path. If the container runs as root
    we change ownership of files written to the host to this
    user and group once everything is finished.

    input: path to get owner and group for
    output: user name, user id, group name, group id
    """
#    login_user_name = os.getlogin()
#    print("login user name:", login_user_name)
    stat_info = os.stat(path_to_check)
    uid = stat_info.st_uid
    gid = stat_info.st_gid
    print("user id:",uid, "group id:", gid)
 
#    user_name = getpass.getuser()
#    print("user name:", user_name)
#    group_id = os.getegid()
#    print("group id1:", group_id)
    user_name = pwd.getpwuid(uid).pw_name
    group_name = grp.getgrgid(gid).gr_name
    print("user name:", user_name)
    print("group name:", group_name)
    return user_name, uid, group_name, gid

def __main__(args):
    """
    """
    start_time = time.time()

    options = parse_arguments()

    #get the user and group information for a path on the host. If the container
    #runs as root then all files written onto the host file system will be owned 
    #by root. To avoid this we chown any written files to the user.
    #We assume the user is the owner of the ouput directory in this case.
#    user_name, user_id, group_name, group_id = get_user_group_for_a_path(options.output_dir)

    yaml_substitute_values = collections.defaultdict(str)

    bed_file_str = "".join(options.bed_file) 
    yaml_substitute_values['bed_file'] = bed_file_str

    output_dir_str = "".join(options.output_dir)
    yaml_substitute_values['output_dir'] = output_dir_str

#    if options.sample_files:
#        sample_file_str = ",".join(options.sample_files)
#        yaml_substitute_values['sample_files'] = sample_file_str

    if options.normal_sample_files:
#        print("normal files:", options.normal_sample_files)
        normal_file_str = ",".join(options.normal_sample_files)
#        print("items", normal_file_list)
        yaml_substitute_values['normal_sample_files'] = normal_file_str

    if options.tumor_sample_files:
        tumor_file_str = ",".join(options.tumor_sample_files)
        yaml_substitute_values['tumor_sample_files'] = tumor_file_str

#    for key, item in yaml_substitute_values.iteritems():
#        print("yaml key and item:\n",key, item)

    # initialize structural variant substitution to nothing
    # in case structural variant calling is not requested
    print("workflow:", options.workflow)

#    if 'germline_variant_calling' in options.workflow:
#        workflow_template = get_germline_variant_template()
    
    if 'cancer_variant_calling' in options.workflow:
        workflow_template = get_cancer_variant_template()
    
    yaml_substitute_values['svcaller_info'] = ""
    if 'structural_variant_calling' in options.workflow:
        workflow_template = get_cancer_variant_template()
        yaml_substitute_values['svcaller_info'] = 'svcaller: [cnvkit, lumpy, delly]' 

    workflow_to_run = ""   
    try:
        workflow_to_run = workflow_template.substitute(yaml_substitute_values)
    except KeyError, err:
        print('ERROR:', str(err), file=sys.stderr) 



    print("workflow template:\n",workflow_to_run)

    if options.GATK_file:
        # if the gatk dir /tmp/gatk exists and is not empty then the GATK may already
        # be installed so do not reinstall the GATK 
        if (os.path.isdir('/tmp/gatk/') and (len(os.listdir('/tmp/gatk/')) != 0)):
            print("WARNING: The container's GATK directory is not empty,"
                "skipping GATK install!", file=sys.stderr)
        else:
            cmd = "gatk-register " + options.GATK_file
            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            except subprocess.CalledProcessError as err:
                print(err.output, file=sys.stderr)
            else:
                print("GATK file ", options.GATK_file, " registered")

    #install the genome data if the provided data directory is empty
    #TODO: check if particular genomes need to be installed...allow
    #user to specify additional genomes to install.
    #The data directory should be mounted to the /mnt/biodata/
    #directory inside the container, so we can look at /mnt/biodata/
    #to see if it exists or is empty.
    if os.path.exists('/mnt/biodata/'):
        # if the data dir is empty then we can safely download
        # genome data
        if len(os.listdir('/mnt/biodata/')) == 0:
            #create the subdirectory which will hold the location
            #files and bcbio system YAML
            os.makedirs('/mnt/biodata/galaxy')
            #change the ownership of the directory to the  user so when the bcbio_nextgen 
            #command is executed as the user below the files in it can be read or edited
#            os.chmod('/mnt/biodata/galaxy', 0666)
#            os.chown('/mnt/biodata/galaxy', user_id, group_id)

            #create the subdirectory that will hold the reference genomes data
            os.makedirs('/mnt/biodata/genomes')
            #change the ownership of the directory to the  user so when the bcbio_nextgen 
            #command is executed as the user below the files in it can be read or edited
#            os.chmod('/mnt/biodata/genomes', 0666)
#            os.chown('/mnt/biodata/genomes', user_id, group_id)

            #output the bcbio system YAML file to /mnt/biodata/galaxy/
            #This contains default values that can be edited by the user later 
            #after the current run to reflect the resources on the system. 
            sys_yaml_substitute_values = collections.defaultdict(str)
            #TODO: add correct system mem to substitute values?
            #set the default number of cores in the system
            sys_yaml_substitute_values['core_count'] = str(16)
            #try to get the true number of cores in the system
            #and insert that into the bcbio system YAML
            core_count = 16
            try:
                core_count = multiprocessing.cpu_count()
            except: 
                print('WARNING: Could not get actual core count; using default:', core_count, file=sys.stderr) 

            sys_yaml_substitute_values['core_count'] = str(core_count)
            bcbio_system_yaml_template = get_bcbio_system_template()
            bcbio_system_yaml = ""  
            #replace the placeholders in the system yaml template 
            try:
                bcbio_system_yaml = bcbio_system_yaml_template.substitute(
                                                   sys_yaml_substitute_values)
            except KeyError, err:
                print('ERROR:', str(err), file=sys.stderr) 

            print("bcbio system YAML:", bcbio_system_yaml)

            with open("/mnt/biodata/galaxy/bcbio_system.yaml","w+") as bcbio_system_file:
                bcbio_system_file.write(bcbio_system_yaml)

            #change the ownership of the file to the  user so when the bcbio_nextgen 
            #command is executed as the user below the files in it can be read or edited
#            os.chown(bcbio_system_file.name, user_id, group_id)
            #os.chmod(bcbio_system_file.name, 0666)

            #run the command to download the data for the reference genome and
            #location files. The reference genome will be put in /mnt/biodata/genomes/
            #and the location files in /mnt/biodata/galaxy/
            #These directories should point to directories on the host that have lots
            #of space and setup by mounting a volume to /mnt/biodata in the 'docker run'
            # command. I.e. docker run -v $(pwd)/data/:/mnt/biodata/ ...
            cmd = ["bcbio_nextgen.py", 'upgrade', '--data', '--genomes', 'GRCh37', '--aligners', 'bwa']


            #run the bcbio_nextgen.py command using the bcbio script that runs the bcbio-nexgen.py
            #program as the user we specify so that ouput files are owned by that user and not the 
            #user the container is running as which could a different user
#            cmd = ["/sbin/createsetuser", user_name, str(user_id), group_name, str(group_id),
#            'bcbio_nextgen.py', 'upgrade', '--data', '--genomes', 'GRCh37', '--aligners', 'bwa'] 

            print("command to run:\n",cmd)
            output = subprocess.call(cmd)
            print("data download output is:\n", output)

        else:
            print("WARNING: The data directory is not empty, skipping data download!", 
                     file=sys.stderr)
    else:
        print("ERROR: The provided data directory is not a directory, skipping"
               "data download!", file=sys.stderr)
 

    #create a temporary file in which to store the YAML template so the bcbio-nextgen
    #script can read it. The file will be deleted after the command is run. 
#     with tempfile.NamedTemporaryFile(dir=options.work_dir, delete=False) as workflow_yaml_file:
    with tempfile.NamedTemporaryFile(delete=False) as workflow_yaml_file:
        workflow_yaml_file.write(workflow_to_run)    
    #os.chmod(workflow_yaml_file.name, 0644)
    #change the ownership of the file to the  user so when the bcbio_nextgen 
    #command is executed as the user below the files in it can be read or edited
#    os.chown(workflow_yaml_file.name, user_id, group_id)


    #run the workflow as the user and group that owns the output directory
    #run the bcbio_nextgen.py command using the bcbio script that runs the bcbio-nexgen.py
    #program as the user we specify so that ouput files are owned by that user and not the 
    #user the container is running as which could a different user
#    cmd = ["/sbin/createsetuser", user_name, str(user_id), group_name, str(group_id), 
#           'bcbio_nextgen.py', workflow_yaml_file.name, "-n", str(options.num_cores)]
    #run the workflow
    cmd = ["bcbio_nextgen.py", workflow_yaml_file.name , "-n", str(options.num_cores)]
     
    print("command to run:\n",cmd)
    output = subprocess.call(cmd)
    print("workflow output is:\n", output)

    #delete the temporary file used to hold the YAML template
    workflow_yaml_file.close()

    #TODO??? change the owner and group of files written to the host to
    #the owner and group captured at the beginning of the script
    #we assume these are the files in the working and output directories
        
    

    print("----- %s seconds -----" % (time.time() - start_time), file=sys.stderr)

if __name__=="__main__":
     sys.exit(__main__(sys.argv))

