FROM bcbio/bcbio:latest

MAINTAINER Walt Shands jshands@ucsc.edu

WORKDIR ./work

USER root

COPY UCSCbcbioTool.py /usr/local/bin
RUN chmod a+x /usr/local/bin/UCSCbcbioTool.py

# switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu

# bcbio_nextgen creates tmp directories in /home/ubuntu so we need to create 
# /home/ubuntu in the root as owned by ubuntu before switching to user ubuntu
RUN mkdir /home/ubuntu
RUN chown ubuntu:ubuntu /home/ubuntu

RUN chmod a+rw /mnt/biodata

USER ubuntu

#Use VOLUME to expose the /mnt/biodata/ directory on the host. This is done to make
#the the /mnt/biodata/ directory writable from inside the container, when the
#container is launched by cwltool. By default the file system in the container
#is set to readonly by cwltool. Bcbio requires that symlinks /mnt/biodata/galaxy
#and /mnt/biodata/genomes point to reference data on the host. This data is 
#typically downloaded by bcbio to a directory with subdirectories 'galaxy' and
#'genomes'. The script inside of this container will create the symlinks
#but it needs /mnt/biodata/ to be writable.

#This data directory will default to the current working directory on the host
#when cwltool does 'docker run'; cwltool sets the current working directory
#to a subdirectory under /tmp. The user can control this by setting TMPDIR on 
#the host environment. The user should definately set this environment variable 
#since /tmp may not have enough space for the reference and intermediate files.

#If the container is run by the user directly exexcuting 'docker run...', and
#the user passes in a data directory where the reference data is already located,
#or an empty directory, in which case the script in will call bcbio to download
#the reference data to that directory, the symlinks will be created.
#Alternatively if the user passes in a path and name of the
#tarred up reference data directory, in which case the data directory will be 
#untarred in the current working directory, the symlinks will be setup.

#NOTE: WHEN RUNNING WITH DOCKSTORE OR CWLTOOL IT IS HIGHLY SUGGESTED THAT THE 
#USER SET THE HOST ENVIRONMENT VARIABLE 'TMPDIR' TO A DIRECTORY WITH ENOUGH 
#SPACE TO ACCOMODATE VERY LARGE FILES; E.G. THE REFERENCE DATA DOWLOADED BY 
#BCBIO CAN CONSUME > 15G
 
VOLUME /mnt/biodata/


