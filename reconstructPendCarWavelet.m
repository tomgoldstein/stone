%  This script creates a stream of Stone measurements from the PendCar
%  high-speed camera dataset, and then reconstructs it.

%%  Parameters that control reconstruction : chosen by user
Nr = 256;  % DMD resolution
dataPerFrame = 4096; %round(Nr*Nr/20);  % 2048  How many samples to use in a frame - i.e. the width of the data window
shiftPerFrame = 4096;% 256/512 round(dataPerFrame/16); % /16for128 % How many samples does the data window shift between frames
mu = .01; % 150,  The weight of the L1 term in the variational reconstruction  
Nf = 100;  %  MAXIMUM number of frames to reconstruct

regularizer = 'wavelets2'; 
%regularizer = 'wavelets3'; 
%regularizer = 'dct3'; 

solver = 'fasta'; 
%solver = 'cosamp'; 


order = createOrderingData(Nr,'full'); % choose either 'full' random, or 'semi' random

%%  Parameters that control data sampling : chosen by user
numImages = 256;   % number of images from PendCar to read from
samplesPerImage = 256;  % stone samples to pull from each image
[ data, finalImage ] = createMeasurementsFromImages( 'PendCar_lowres', ...
                    numImages, samplesPerImage, order);



%% Parameters that we can calculate
% Number of data
Nd = size(data,1);  % number of data
% Number of frames
Nf = min(Nf,floor((Nd-dataPerFrame)/shiftPerFrame)+1);
% Number of pixels
Np = Nr*Nr;  

%  The time dimension must have power of two length so that boundary
%  reflections don't make the transform non-orthogonal.
if strcmp(regularizer,'wavelets3')
    Nf = 2^floor(log2(Nf));
end
    

fprintf('Reconstructing:\n\tResolution = %d\n\tFrames = %d\n\tSamples = %d (%3.2d)%%\n\tShift = %d\n',Nr,Nf,dataPerFrame,round(dataPerFrame/Np*100000)/1000,shiftPerFrame);

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
% %%  Subtract means to convert the data from 0/1 to +1/-1
% for f=1:Nf
%     avZ = sum(b(:,f))/sum(R(:,f));
%     avSTO = avZ/(Nr+1);
%     b(:,f) = b(:,f) + (avSTO-avZ); 
% end
fprintf('%f secs\n',toc);

%  Options for FASTA
opts = [];
opts.verbose = 2;
%opts.accelerate = true;

if strcmp(regularizer,'wavelets2')
    level = floor(log2(Nr));
    [C,S] = wavedec2(zeros(Nr,Nr),level, 'haar');  %  setup wavelet transform
    %  The measurement operator
    M = @(haarVecs) R.*STO(haar2ToImage( haarVecs, S, order ));
    Mt = @(stoneVecs) imageToHaar2( STO(R.*stoneVecs), level, order );
    x0 = zeros(Np,Nf);
    opts.tol = 1e-4;
elseif strcmp(regularizer,'wavelets3')
    level = min(floor(log2(Nr)),floor(log2(Nf)));
    [S] = wavedec3(zeros(Nr,Nr,Nf),level, 'haar');
    %  The measurement operator
    M = @(haarVecs) R.*STO(haar3ToImage( haarVecs, S, order ));
    Mt = @(stoneVecs) imageToHaar3( STO(R.*stoneVecs), level, order );
    x0 = imageToHaar3(zeros(Np,Nf),level,order);
    opts.tol = 1e-4;
elseif strcmp(regularizer,'dct3')
    %  The measurement operator
    M = @(d) R.*STO(dct3ToImage( d, order ));
    Mt = @(stoneVecs) imageToDct3( STO(R.*stoneVecs), order );
    x0 = imageToDct3(zeros(Np,Nf),order);
    opts.tol = 1e-4;
else
    assert(false,'invalid regularizer');
end

if strcmp(solver,'cosamp')
    vec = @(x) x(:);
    xshape = @(x) reshape(x,size(x0));
    dshape = @(x) reshape(x,size(b));
    M = @(x) vec(M(xshape(x)));
    Mt = @(x) vec(Mt(dshape(x)));
    K = sum(R(:));
    b = b(:);
end


%  SOLVE
fprintf('Reconstructing...\n'); tic;
% setup fasta
f = @(x) 0.5*norm(x-b,'fro')^2;
gradf = @(x) x-b;
g = @(x) mu*norm(x(:),1);
proxg = @(x, tau)  sign(x).*max(abs(x)-mu*tau,0);

%% Call the solver to minimize L1 and recover image from sub-sampled data
if strcmp(solver,'fasta')
    tic;
    [x,outs ] = fasta(M, Mt, f, gradf, g, proxg, x0, opts );
    time = toc;
    fprintf('FASTA took %f seconds\n',time);
elseif strcmp(solver,'cosamp')
    tic;
    [x,r,normR,residHist, errHist] = CoSaMP( {M,Mt}, b, K, [], opts );
    x = xshape(x);
    time = toc;
    fprintf('CoSaMP took %f seconds\n',time);
else
    assert('invalid solver option');
end

if strcmp(regularizer,'wavelets2')
    recon = haar2ToImage( x, S, order );
elseif strcmp(regularizer,'wavelets3')
    recon = haar3ToImage( x, S, order );
elseif strcmp(regularizer,'dct3')
    recon = dct3ToImage( x, order );
else
    assert(false,'invalid regularizer');
end

%disp('Refining Results');
%[u,y,tau,rp,rd,A,At] = pdhg_video_L0( R, b, mu, order,u,y);

%% Print how many frames
fprintf('Reconstructed %d frames\n', Nf);

%% Create an array of 2d frames from the column vectors by re-shaping them
frames = zeros(Nr,Nr,Nf);
for f = 1:Nf
    frames(:,:,f) = nestedVectorToImage(recon(:,f),order);
end

frames = cleanImage(frames);

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


filename = ['DMD_reconstructions/pendcar/waveletfbf_Nr_',num2str(Nr) ,'_Nd_' ,num2str(length(data)) ,'_Ndf_', num2str(dataPerFrame),'_shift_',num2str(shiftPerFrame),'_Nf_',num2str(Nf), '_mu_',num2str(mu) ];
imdata = permute(frames,[1 2 4 3]);
imwrite(imdata,[filename '.gif'],'DelayTime',0,'LoopCount',inf);

imwrite(frames(:,:,end/2)/256,['DMD_reconstructions/pendcar/wavelet_sample_frame.png'],'mode','lossless');

save([filename '.mat'],'frames');
