
# r442_tidyverse_seurat5

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Rstudio
 - R 4.4.2
 - Misc packages for data manipulation, figures, reports
 - Seurat 5

# Note : It is a good practice to put the name of the group and project into the image name while compiling
# in order to be able to associate easily the group and project name to the image

## Build (in directory where the Docker file is)

docker build . -t evlab_labeltransfer_r442_seurat5

## Save
## Save the image to the project folder into the experiment 02_Container folder (in lab side)

docker save evlab_labeltransfer_r442_seurat5 | gzip > /mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/labelTransferBenchmark/01_SeuratEvaluation/02_Container/evlab_labeltransfer_r442_seurat5.tar.gz

## Run Rstudio
## Command to start a container based on the image.
## You have to change "yourPass" to a simple password (no space, no special characters).
## You can also change the mapped port number (here 9090).
docker 
Give user details to internal script which sets user and permissions:

```
docker run -d --name Rstud -p 9090:8787 -e PASSWORD=pw -e DEFAULT_USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt evlab_labeltransfer_r442_seurat5
```

Then connect to the machine running docker (localhost) on mapped port (9090):
http://127.0.0.1:9090



