# inference-wound-healing
This repository provides the code and processed data used in the manuscript: "Inference of weak-form partial differential equations describing migration and proliferation mechanisms in wound healing experiments on cancer cells."

## Directory structure

The code and results for Sec 4.1 and sec 4.2 are located in directories: "wound_healing_jin/" and "wound_healing_tram/", respectively. 
Each directory consists of 4 subdirectories: 
1. data/: directory contains the smooth density data used for inference. 
2. FenicsCode/: directory contains the code used for data-driven inference. All the codes are based on the Fenics library (see installation instructions here https://fenicsproject.org/download/archive/)
The codes are to be run sequentially from Step1-5. 
3. results/: directory contains the results generated while running the scripts in FenicsCode
4. postprocessing/: directory contains the MATLAB plotting scripts used to generate the paper figures. Use codes labeled plot_* to plot the manuscript figures.

wound_healing_jin/jinModelFenics contains the parameters extracted from Jin et.al. 

## Docker usage instructions for fenics code. 

Step 1: Install Docker

Step 2: Fix the image to be pulled

export image_name=quay.io/fenicsproject/dev-env
export tag=latest
Step 3: Pull the image to your local system

docker pull  $image_name:$tag
Step 4: Run a container named fenics from the image on your local system

docker run --name fenics -dit -w /home/fenics/shared -v $(pwd):/home/fenics/shared $image_name:$tag
Step 5: Go into the container

docker exec -it fenics bash
Coming out of container: ctrl+D

Step 6: Run the relevant python file.

Other useful commands:

Removing an image (You will not need this unless you used the wrong image)
docker rmi $image_name:$tag
Removing a container named fenics
docker rm fenics
Checking status of all containers
docker ps -a
Checking status of all images
docker ps -i
Running with privileges to get mount on
docker run --name fenics -dit --privileged -w /home/fenics/shared -v $(pwd):/home/fenics/shared $image_name:$tag

