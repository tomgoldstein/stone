 %  This script loads different DMD datasets, and calls the reconstruction
%  method.  Results are saved in the 'DMD_reconstructions' folder.


%% Reconstruct 'car' datasets
name = 'cars2';
type = 'l2';
load 'dmd_measurements/car/128/data2.mat';

Nr = 128;  % DMD resolution
dataPerFrame = 1024;
reconstructMovingDMDMeasurements; 

load 'dmd_measurements/car/256/32x/data1.mat';
Nr = 256;  % DMD resolution
dataPerFrame = 1024;
reconstructMovingDMDMeasurements; 


%% Reconstruct 'hand' datasets

name = 'hand2';
type = 'l2';

load 'dmd_measurements/hand/128/data2.mat';
Nr = 128;  % DMD resolution
dataPerFrame = 1024;
reconstructMovingDMDMeasurements; 

load 'dmd_measurements/hand/256/data2.mat';
Nr = 256;  % DMD resolution
dataPerFrame = 1024;
reconstructMovingDMDMeasurements; 
