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
    The sample fastq or BAM files: These should be separated by a space. 
        Required for germline variant calling.
    The tumor fastq or BAM files: These should be separated by a space. 
        Required for cancer variant calling.
    The normal fastq or BAM files: These should be separated by a space.
        Required for cancer variant calling.
    Total available cores. This tells bcbio how many total cores to use. 
               http://bcbio-nextgen.readthedocs.io/en/latest/contents/parallel.html
    The path and file name of the GATK tools. If provided the GATK tools will be installed
          in the container for use in the workflow. (required????) 
    The name of the workflow to run. If not provided cancer variant calling only is run.
    The path to and name of the BED file. Required.
    The path to the tar file of the bcbio reference genome directory or a path
        to where the bcbio reference genome directory already exists or where the 
        reference genomes should be download. 
        For Dockstore this must be the
        tar file of the reference genome directory or nothing, in which case
        the reference genomes are downloaded to the current working directory 
        specified by the users TMPDIR environment variable, or if that is not
        set then then cwltool usually sets the cwd to a subdirectory of /tmp
    The path of the output directory. Default is cwd/final. 

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

#    worflow_input_switches = parser.add_mutually_exclusive_group()
    germline_input_switches= parser.add_argument_group('germline input switch')
    germline_input_switches.add_argument('-s', '--sample_files', type=str, action='append',
                         help="Input file(s) for processing. Multiple path/file"  
                         "names can be provided separated by spaces. These could" 
                         "correspond fastq files for paired end reads." )

    cancer_input_switches = parser.add_argument_group('cancer input switches')
    cancer_input_switches.add_argument('-t', '--tumor_sample_files', type=str,
                         action='append',
                         help="Input file(s) for processing. Multiple path/file"
                         "names can be provided separated by spaces. These could"
                         "correspond fastq files for tumor paired end reads." )
    cancer_input_switches.add_argument('-n', '--normal_sample_files', type=str,
                         action='append', help="Input file(s) for processing."
                         "Multiple path/file names can be provided separated by"
                         "spaces. These could correspond fastq files for normal"
                         "paired end reads." )

    parser.add_argument('-c', '--num_cores',  type=int, default=16,  
                         help="number of cores to use for processing.") 

    parser.add_argument('-g', '--GATK_file', type=str, help='Path to GATK file,' \
                        'e.g. /path/to/GenomeAnalysisTK.tar.bz2.')
 
    parser.add_argument('-W','--workflow', choices=['cancer-variant-calling',
                         'germline-variant-calling','structural-variant-calling'], 
#                         default=['cancer-variant-calling'], 
                         action='append',
                         required=True, 
                         help="The name of the workflow to run. Default is cancer" \
                         " variant calling.")

    parser.add_argument('-b', '--bed_file', type=str, required=True, 
                          help='Input BED file.')

    reference_genomes_group  = parser.add_mutually_exclusive_group()
    reference_genomes_group.add_argument('-d', '--data_dir', type=str, 
                       help='Directory where reference genome files are or should'
                       'be downloaded.' )
    reference_genomes_group.add_argument('-f', '--data_file', type=str, help='Path to reference genomes tar file.') 
 
    parser.add_argument('-o', '--output_dir', type=str, 
                       help='Directory where output files should be written. Default is <cwd>/final' )

    options = parser.parse_args()


    print("workflows:", options.workflow)
    if ('germline-variant-calling' not in options.workflow  and 
       'cancer-variant-calling' not in options.workflow and 
       'structural-variant-calling' in options.workflow):
        parser.error('Structural variant calling must be run with germline or cancer variant calling')

    if ('germline-variant-calling' in options.workflow  and 
       'cancer-variant-calling' in options.workflow):
        parser.error('Cancer variant calling cannot be run with germline variant calling')

    if 'germline-variant-calling' in options.workflow  and options.sample_files is None:
        parser.error('Sample files switch must be used for germline-variant-calling')

    if ('cancer-variant-calling' in  options.workflow  and 
                 ((options.normal_sample_files is None) or (options.tumor_sample_files is None))):
        parser.error('Normal and tumor input file switches must be used for cancer-variant-calling')

    if(((options.normal_sample_files is not None) or (options.tumor_sample_files is not None)) and
        options.sample_files is not None):
         parser.error('Normal and tumor input file switches cannot be used with sample input switch')

    return (options)


def get_bcbio_system_template():
    """
    This templated describes the system resources to the bcbio Docker container
    so that it can allocate the correct resources to each thread.
    """
    bcbio_system_template = string.Template("""
#
resources:
#  tmp:
#    dir: $temp_dir
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
#  galaxy_config: $galaxy_config_path
    """)

    return bcbio_system_template



def get_germline_variant_template():
    """
    This is the template that describes to bcbio the tools to used in the cancer
    variant calling workflow and if structural variant calling is also requested
    which tools to use. The input tumor and normal path and file names are 
    inserted into the template as is the path and file name of the BED file,
    the structural variant calling tools to use if structural variant calling is
    requested and the output directory.
    """
    germline_variant_template = string.Template("""
# See the bcbio-nextgen documentation for a description of the pipeline this was
# derived from 
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines
---
resources:
  tmp:
    dir: $working_dir
    #dir: ./work
upload:
  dir: $output_dir
details:
#   - files: [../input/NA12878_1.fastq.gz, ../input/NA12878_2.fastq.gz]
  - files: [ $sample_files ]
    description: NA12878
    metadata:
      batch: ceph
      sex: female
    analysis: variant2
    genome_build: GRCh37
    algorithm:
      aligner: bwa
      align_split_size: 5000000
      mark_duplicates: true
      recalibrate: false
      realign: false
      variantcaller: [freebayes, gatk-haplotype, platypus, samtools]
      remove_lcr: true
#      validate: ../input/GiaB_v2_19.vcf.gz
#      validate_regions: ../input/GiaB_v2_19_regions.bed
    """)
    return germline_variant_template



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
resources:
  tmp:
    #dir: ./work 
    dir: $working_dir
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

    #get the current working directory. If the container is 
    #run by cwltool alone or by cwltool via a Dockstore 
    #run command then all of the intermediate and results files will be written
    #to the current working directory. It is important for the user in this case
    #to set the host environment variable TMPDIR to a directory will plenty of 
    #space in which to write, otherwise cwltools will create a temporary directory,
    #most likely under /tmp, for this purpose, and there may not be enough disk
    #space for the pipelines to use.
    cwd = os.getcwd()
    if not os.path.exists(cwd):
        print("ERROR: The current working directory " + cwd + " cannot be reached," 
              "did you mount a volume to the directory using '-v' on the 'docker run'"
              "command line?", file=sys.stderr)
 
    #if the user has provided a directory where the bcbio reference data has
    #been already downloaded or where it should be downloaded
    if options.data_dir:
        datadir = options.data_dir    
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        #if the data dir does not end in a slash add one
        datadir = os.path.join(datadir, '')
    #if the user has  provided a tar file of the tarred up bcbio reference
    #data 
    elif options.data_file:
        #create a subdirectory in which to extract bcbio reference data
        #This should be the current working directory as set by the user
        #or cwltool as we know it will be writable and we hope that the 
        #user has chosen a directory with enough space
        #TODO??? This should be a unique name???? in case it already exists!!! 
        datadir = cwd + '/data/'
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        #open the bcbio reference data tar in the current working directory
        #Dockstore (and cwltool?) seems to require this as the directory
        #where the input tar file is stored (and all other inputs) are
        #mounted as read only 
        cmd = ['tar','-xzvf', options.data_file, '-C', datadir]
        print("cmd is:", cmd) 
        print("Untarring bcbio data file ", options.data_file, " in ", datadir)
        try:
            subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as err:
            print(err.output, file=sys.stderr)
        else:
            print("Bcbio data file ", options.data_file, " untarred in ", datadir)
    else:
        #create a subdirectory in which to store bcbio reference data. The user 
        #(or cwltool) has not provided a directory for the reference data and no
        #tar file of the reference data has been provided. In this case the reference 
        #has to be downloaded, and the default is to put it in a subfolder under cwd 
        #TODO??? This should be a unique name???? in case it already exists!!! 
        datadir = cwd + '/data/'
        if not os.path.exists(datadir):
            os.makedirs(datadir)
 
      
    print("bcbio reference data dir:", datadir)

    #create a directory for the intermediate files created during the pipeline
    #run. The default is to put them under the current working directory. 
    #This path is substituted into the bcbio workflow YAML so bcbio will write
    #intermediate files to this directory.
    working_dir = cwd + '/work/'
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)    

#    sys_yaml_substitute_values['working_dir'] = working_dir
 
    #get the user and group information for a path on the host. If the container
    #runs as root then all files written onto the host file system will be owned 
    #by root. To avoid this we chown any written files to the user.
    #We assume the user is the owner of the ouput directory in this case.
#    user_name, user_id, group_name, group_id = get_user_group_for_a_path(options.output_dir)

    yaml_substitute_values = collections.defaultdict(str)

    yaml_substitute_values['working_dir'] = working_dir

    bed_file_str = "".join(options.bed_file) 
    yaml_substitute_values['bed_file'] = bed_file_str

    
    #create a directory for the result files created during the pipeline
    #run. The default is to put them under the current working directory. 
    #This path is substituted into the bcbio workflow YAML so bcbio will write
    #intermediate files to this directory.
    output_dir_str = './final'
    if options.output_dir: 
        output_dir_str = "".join(options.output_dir)
    if not os.path.exists(output_dir_str):
      os.makedirs(output_dir_str)    
    yaml_substitute_values['output_dir'] = output_dir_str

    if options.sample_files:
        sample_file_str = ",".join(options.sample_files)
        yaml_substitute_values['sample_files'] = sample_file_str

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

    print("workflow:", options.workflow)

    if 'germline-variant-calling' in options.workflow:
        workflow_template = get_germline_variant_template()
    
    if 'cancer-variant-calling' in options.workflow:
        workflow_template = get_cancer_variant_template()
    
    # initialize structural variant substitution to nothing
    # in case structural variant calling is not requested
    yaml_substitute_values['svcaller_info'] = ""
    if 'structural-variant-calling' in options.workflow:
        yaml_substitute_values['svcaller_info'] = 'svcaller: [cnvkit, lumpy, delly]' 

    #place the input files and intermediate and results files locations, etc.
    #into the bcbio workflow YAML and save the YAML to pass to bcbio
    workflow_to_run = ""   
    try:
        workflow_to_run = workflow_template.substitute(yaml_substitute_values)
    except KeyError, err:
        print('ERROR:', str(err), file=sys.stderr) 

    print("workflow template:\n",workflow_to_run)

    #install the GATK inside the container so bcbio can use the GATK tools
    #The default location to install the tools inside the container is /tmp
    #so check first to see if somehow they are already installed
    if options.GATK_file:
        # if the gatk dir /tmp/gatk exists and is not empty then the GATK may already
        # be installed so do not reinstall the GATK 
        if (os.path.isdir('/tmp/gatk/') and (len(os.listdir('/tmp/gatk/')) != 0)):
            print("WARNING: The container's GATK directory is not empty,"
                "skipping GATK install!", file=sys.stderr)
        else:
            cmd = ["gatk-register", options.GATK_file]
            try:
                #TODO try formatting the command as one string and removing shell=True
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as err:
                print(err.output, file=sys.stderr)
            else:
                print("GATK file ", options.GATK_file, " registered")

    #install the genome data if the provided data directory is empty
    #TODO: check if particular genomes need to be installed...allow
    #user to specify additional genomes to install.

    #The reference data subdirectory 'genomes' must be pointed to by the 
    #/mnt/biodata/genomes directory inside the container, and the same is true 
    #for the subdirectory 'galaxy'. This is a bcbio requirement and is accomplished 
    #by creating symlinks. 

    #In the situation where cwltool runs the container the file system is read
    #only and the the script would not be able to create the symlinks. We avoid
    #this by using the VOLUME /mnt/biodata statment in the Dockerfile and setting 
    #the directory permissions so that /mnt/biodata is writeable by any user.

    #if the call to docker run does not include a -v mount point for the
    #genome reference data on the host then the container will fail to see
    #the directory so keep this test to check for this
    if os.path.exists(datadir):
        # if the data dir is empty then we can safely download genome data
        #if it is not empty we assume the data has already been downloaded
        #and bcbio can use it
        if len(os.listdir(datadir)) == 0:
            #create the subdirectory which will hold the location
            #files and bcbio system YAML
            os.makedirs(datadir + 'galaxy')
            #create the necessary symlink so bcbio can find the galaxy folder on the host
            os.symlink(datadir + 'galaxy', '/mnt/biodata/galaxy')

            #change the ownership of the directory to the  user so when the bcbio_nextgen 
            #command is executed as the user below the files in it can be read or edited
#            os.chmod('/mnt/biodata/galaxy', 0666)
#            os.chown('/mnt/biodata/galaxy', user_id, group_id)

            #create the subdirectory that will hold the reference genomes data
            os.makedirs(datadir + 'genomes')
            #create the necessary symlink so bcbio can find the genomes folder on the host
            os.symlink(datadir + 'genomes', '/mnt/biodata/genomes')

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

            with open('/mnt/biodata/galaxy/' \
                      + "bcbio_system.yaml","w+") as bcbio_system_file:
                bcbio_system_file.write(bcbio_system_yaml)

            #change the ownership of the file to the  user so when the bcbio_nextgen 
            #command is executed as the user below the files in it can be read or edited
#            os.chown(bcbio_system_file.name, user_id, group_id)
            #os.chmod(bcbio_system_file.name, 0666)

            #run the command to download the data for the reference genome and
            #location files. The reference genome will be put in /mnt/biodata/genomes/
            #and the location files in /mnt/biodata/galaxy/
            #These directories should point to directories on the host that have lots
            #of space and this is accomplished via the symlinks created above
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
            #TODO check that the genomes and galaxy folders exist; they should if 
            #datadir points to a bcbio created data directory
            #create the necessary symlink so bcbio can find the genomes folder on the host
            print("setting symlink /mnt/biodata/galaxy to point to:", datadir+'galaxy')
            os.symlink(datadir + 'genomes', '/mnt/biodata/genomes')
            #create the necessary symlink so bcbio can find the galaxy folder on the host
            os.symlink(datadir + 'galaxy', '/mnt/biodata/galaxy')
            print("WARNING: The data directory " + datadir + "is not empty, skipping data download!", 
                     file=sys.stderr)
    else:
        print("ERROR: The provided data directory " + datadir + " cannot be reached," 
              "did you mount a volume to the data directory using '-v' on the 'docker run'"
              "command line?", file=sys.stderr)
 
    #write the project YAML to the current working directory so the user can
    #see exactly what bcbio is running and so bcbio can read it and run it
    with open("bcbio_project.yaml","w+") as bcbio_project_file:
                bcbio_project_file.write(workflow_to_run)
    bcbio_project_file.close()


    #os.chmod(workflow_yaml_file.name, 0644)
    #change the ownership of the file to the  user so when the bcbio_nextgen 
    #command is executed as the user below the files in it can be read or edited
#    os.chown(workflow_yaml_file.name, user_id, group_id)

    #run the workflow
    cmd = ["bcbio_nextgen.py", bcbio_project_file.name , "-n", str(options.num_cores)]
     
    print("command to run:\n",cmd)
    output = subprocess.call(cmd)
    print("workflow output is:\n", output)

    #TODO??? change the owner and group of files written to the host to
    #the owner and group captured at the beginning of the script
    #we assume these are the files in the working and output directories
        
    

    print("----- %s seconds -----" % (time.time() - start_time), file=sys.stderr)

if __name__=="__main__":
     sys.exit(__main__(sys.argv))

