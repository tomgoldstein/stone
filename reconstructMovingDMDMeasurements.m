%   This script creates a video from DMD data
%  You must load a data file from the folder 'dmd_measurements' before you 
%  run this script.  See the script 'reconstructDMDVideos.m', which loads 
%  a variety of datasets, and then calls this script for reconstruction.

%name = 'cars';
%load 'dmd_measurements/car/128/data1.mat';


%%  Parameters to be chosen by the user
%Nr = 128;  % DMD resolution
%dataPerFrame = 2048;  % 1024  How many samples to use in a frame - i.e. the width of the data window
shiftPerFrame = 512;% 128/512 round(dataPerFrame/16); % /16for128 % How many samples does the data window shift between frames
mu = 150; % 150,  The weight of the data term in the variational reconstruction  
Nf = 1000;  %  MAXIMUM number of frames to do


order = createOrderingData(Nr,'full'); % choose either 'full' random, or 'semi' random

%% Parameters that we can calculate
% Number of data
Nd = size(data,1);  % number of data
% Number of frames
Nf = min(Nf,floor((Nd-dataPerFrame)/shiftPerFrame));
% Number of pixels
Np = Nr*Nr;  

fprintf('Reconstructing:\n\tResolution = %d\n\tFrames = %d\n\tSamples = %d (%g)%%\n\tShift = %d\n',Nr,Nf,dataPerFrame,round(dataPerFrame/Np*100000)/1000,shiftPerFrame);

%% Allocate memory - one column for each frame
R = zeros(Np,Nf);
b = zeros(Np,Nf);
rowsToSample = order.samplingOrder;
%rowsToSample = circshift(rowsToSample,-1);

fprintf('Binning...'); tic;
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
fprintf('Reconstructing...\n'); tic;
if exist('type','var') && strcmp(type,'l1')
    [u,outs] = pdhg_video_l1( R, b, mu, order );
    type = 'l1';
else
    [u,outs] = pdhg_video( R, b, mu, order );
    type = 'l2';
end
fprintf('Reconstruction took %f secs\n',toc);


%disp('Refining Results');
%[u,y,tau,rp,rd,A,At] = pdhg_video_L0( R, b, mu, order,u,y);

%% Print how many frames
fprintf('Reconstructed %d frames\n',Nf);

%% Create an array of 2d frames from the column vectors by re-shaping them
frames = zeros(Nr,Nr,Nf);
for f = 1:Nf
    frames(:,:,f) = nestedVectorToImage(u(:,f),order);
end

frames = cleanImage(frames, true); % clean up the results and apply a median filter

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

if ~exist('reconstructions')
    mkdir('reconstructions')
end

filename = ['reconstructions/' name '_' type '_Nr_',num2str(Nr) ,'_Nd_' ,num2str(length(data)) ,'_Ndf_', num2str(dataPerFrame),'_shift_',num2str(shiftPerFrame),'_Nf_',num2str(Nf), '_mu_',num2str(mu) ];
imdata = permute(frames,[1 2 4 3]);
imwrite(imdata,[filename '.gif'],'DelayTime',0,'LoopCount',inf);

save([filename '.mat'],'frames');
