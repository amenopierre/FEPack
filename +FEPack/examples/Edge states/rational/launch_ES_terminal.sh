#!/bin/bash

matlab -nodisplay -nosplash -r "addpath(genpath('../../FEPack')); addpath(genpath('../+FEPack')); main_edge_states_rational($1, $2); exit"