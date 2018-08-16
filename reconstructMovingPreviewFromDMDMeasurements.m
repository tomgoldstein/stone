%   This script creates a video from DMD data using preview/direct/nyquist
%   reconstructions of the STO transform.  Results will get saved to a
%   folder called 'DMD_reconstructions'.
%
%  You must load a data file from the folder 'dmd_measurements' before you 
%  run this script.

name = 'cars'; %  cars/hand,  Name of the dataset, used to create filenames to save results
%load 'dmd_measurements/hand/256/data1.mat';
load 'dmd_measurements/car/128/data1.mat';


%%  Parameters to be chosen by the user
Nr = 128;  % DMD resolution
Nprev = 32; % preview resolution
shiftPerFrame = 512;% 128/512, How many samples does the data window shift between frames

Nf = 1000;  %  MAXIMUM number of frames to do
mult = 1.0; %  scale up image to make it brighter before saving to gif 
filter = Nprev>32;  %  use a median filter to clean up video?

%% Parameters that are calculated from the options set by the user
dataPerFrame = Nprev*Nprev;  % 1024  How many samples to use in a frame - i.e. the width of the data window
order = createOrderingData(Nr,'full'); % choose either 'full' random, or 'semi' random
% Number of data
Nd = length(data);  % number of data
% Number of frames
Nf = min(Nf,floor((Nd-dataPerFrame)/shiftPerFrame));
% Number of pixels
Np = Nr*Nr;  

fprintf('Reconstructing Previews (%d):\n\tResolution = %f\n\tFrames = %f\n\tSamples = %f (%d)%%\n\tShift = %f\n',Nprev,Nr,Nf,dataPerFrame,round(dataPerFrame/Np*100000)/1000,shiftPerFrame);
%% Allocate memory - one column for each frame
R = zeros(Np,Nf);
b = zeros(Np,Nf);
rowsToSample = order.samplingOrder;

fprintf('Binning Data into Frames...'); tic;
%%  Bin the data into the columns of b.  Every column contains the transform coefficients for a single frame
start = 1;
stop = dataPerFrame;
for f = 1:Nf
    rowsInThisFrame = rowsToSample(1:dataPerFrame);
    R(rowsInThisFrame,f) = 1;
    b(rowsInThisFrame,f) = data(start:stop);
    start = start+shiftPerFrame;
    stop = stop+shiftPerFrame;
    rowsToSample = circshift(rowsToSample,-shiftPerFrame);
end
fprintf('%f secs\n',toc);

fprintf('Converting to +/-...'); tic;
%%  Subtract means to convert the data from 0/1 to +1/-1
for f=1:Nf
    avZ = sum(b(:,f))/sum(R(:,f));
    avSTO = avZ/(Nr+1);
    b(:,f) = b(:,f) + (avSTO-avZ); 
end
fprintf('%f secs\n',toc);



%  SOLVE
frames = zeros(Nprev,Nprev,Nf);
fprintf('Reconstructing Previews...\n'); tic;
for f = 1:Nf
    frames(:,:,f) = makePreview(R(:,f),b(:,f),order);
end

fprintf('Reconstruction took %f secs\n',toc);

fprintf('Median Filtering...\n'); tic;
frames = cleanImage(frames,filter);
frames = flipud3(frames);
fprintf('Median Filtering took %f secs\n',toc);


%%  Display 4 frames
subplot(2,2,1);
imagesc(frames(:,:,1))
subplot(2,2,2);
imagesc(frames(:,:,2))
subplot(2,2,3);
imagesc(frames(:,:,round(Nf/2)))
subplot(2,2,4);
imagesc(frames(:,:,Nf))

colormap gray;

small = min(frames(:));
big = max(frames(:));
frames = (frames-small)*255/(big-small);

%%  Write the results to an animated gif file
filename = ['DMD_reconstructions/' name '_' num2str(Nr) '_' num2str(Nprev) '_Nf_' num2str(Nf) ] ;
imdata = permute(mult*frames,[1 2 4 3]);
imwrite(imdata,[filename '.gif'],'DelayTime',0,'LoopCount',inf);

save([filename '.mat'],'frames');
