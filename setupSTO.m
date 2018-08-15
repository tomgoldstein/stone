% Set the path to include all subfolders, and compile STO_fast.c

addpath('dataSimulation');
addpath('ordering');
addpath('solvers');
addpath('transforms');
addpath('changeDetection');

cd transforms/
mex STO_fast.c
cd ..