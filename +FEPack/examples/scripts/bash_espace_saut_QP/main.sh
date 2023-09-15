#!/bin/bash

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;92m'
NC='\033[0m' # No Color

# Folder for outputs 
chemin='/home/pierre/Documents/Etudes/these/codes/FEPack/+FEPack';

cheminDonnees='/home/pierre/Documents/Etudes/these/codes/FEPack/+FEPack/outputs';
cd $chemin

# Delete outputs/* if requested
while getopts ":r" opt; do
  case $opt in
    r)
      echo -e "${GREEN}-r was triggered: deleting $cheminDonnees/*${NC}" >&2
      rm $cheminDonnees/*
      ;;
    \?)
      echo -e "${RED}Invalid option: -$OPTARG. Exit code${NC}" >&2
      exit 1
      ;;
  esac
done

# Parameters
numFloquetPoints=128;  # Total number of Floquet points 
sizeFloquetClust=32;   # Number of Floquet points per cluster
(( numClusters=numFloquetPoints/sizeFloquetClust )); # Number of Clusters
numbersNodes=(5 10);  # List of numbers of nodes

# Initialization (parallel wrt numNodes)
for numNodes in "${numbersNodes[@]}"
do
  # Generate the meshes (if not done yet) and compute the FE matrices
  matlab -nosplash -nodesktop -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); init_script_transmission($numFloquetPoints, $numNodes, '$cheminDonnees'); quit;" &
done
wait

for numNodes in "${numbersNodes[@]}"
do
  # Solve the waveguide problems
  for (( idCluster=0; idCluster<$numClusters; idCluster++ ))
  do
    # Show progress
    echo -e "${RED}/////////////////  Cluster $idCluster  /////////////////${NC}"

    # Solve the waveguide problems associated to the current cluster
    matlab -nosplash -nodesktop -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); solve_TFB_waveguide($idCluster*$sizeFloquetClust+1, ($idCluster+1)*$sizeFloquetClust, $numNodes, '$cheminDonnees'); quit;" &
  done
  wait # for all the waveguide problems to be solved
  
  # Conclusion
  matlab -nosplash -nodesktop -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); numNodes = $numNodes; cheminDonnees = '$cheminDonnees'; end_script_transmission; quit;"

  # Remove outputs
  rm $chemin/outputs/TFBU_*
done