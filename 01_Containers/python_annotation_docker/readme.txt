
#This image contains:
# Python and other stuff to run MarkerCount


# ######################
     COMPILE THE IMAGE
# ######################

docker build . -t pyth_ann


#######################
     SAVE IMAGE
#######################

docker save pyth_ann | gzip > pyth_ann.tar.gz


# ######################
     RUN THE IMAGE
# ######################

docker run -d --name pyth_ann -p 8989:8989 -w /mnt -v /mnt/DOSI:/mnt/DOSI -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -e PASSWORD=pw pyth_ann


########################
   CREATE SINGULARITY
########################

docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp/test:/output \
--privileged -t --rm \
quay.io/singularity/docker2singularity \
myDockerImageID

# ATTENTION : ImageID est dc85080d4574
# L'image est alors dans : carer@ito:/tmp/test$


mv ID.sif /mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/labelTransferBenchmark/01_SeuratEvaluation/02_Container/02_Singularity/09-07-

# mv source_path destination




