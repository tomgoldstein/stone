%  Test the compressive video reconstruction using synthetic data.
%  Compressive measurements are created using the +/-1 STO transform.  The
%  data is then reconstructed using the pdhg method.

%%  Set up the problem
n = 256; % use n by n image
N = n*n; % number of pixels
Order = createOrderingData(n,'full'); % create the ordering to map into a vector 
samplesPerFrame = round(n*n/20);
numberOfFrames = 5;
mu = 50;


%%  Create synthetic compressive measurements for the frame
%   'Frames' is a cell array. Each entry is a struct with a row-selector
%   and vector of measurements for a frame
[R B] = createCompressiveMeasurements(numberOfFrames, samplesPerFrame, Order);

B = B + .00*randn(size(B));% Add  noise
B = R.*B;    % Zero-out the unknown entries of b (they got filled in by the noise adding step above)

% Call the pdhg solver to reconstruct
[u,outs] = pdhg_video_l1( R, B, mu, Order );


%  Convert vectorized frames into 2D images
frames = zeros(n,n,numberOfFrames);
for f = 1:numberOfFrames
    frames(:,:,f) = nestedVectorToImage(u(:,f),Order);
end

frames = cleanImage(frames);

%%  Display 4 frames
subplot(4,2,1);
imagesc(frames(:,:,1))
subplot(4,2,2);
imagesc(frames(:,:,2))
subplot(4,2,3);
imagesc(frames(:,:,round(numberOfFrames/2)))
subplot(4,2,4);
imagesc(frames(:,:,numberOfFrames))

%  Create lower resolution data by averaging blocks
[r,b, order ] = reduceDataResolution(R,B, Order);
[u,outs] = pdhg_video_l1( r, b, mu/2, order );

frames2 = zeros(n/2,n/2,numberOfFrames);
for f = 1:numberOfFrames
    frames2(:,:,f) = nestedVectorToImage(u(:,f),order);
end

frames2 = cleanImage(frames2);

%%  Display 4 frames
subplot(4,2,5);
imagesc(frames2(:,:,1))
subplot(4,2,6);
imagesc(frames2(:,:,2))
subplot(4,2,7);
imagesc(frames2(:,:,round(numberOfFrames/2)))
subplot(4,2,8);
imagesc(frames2(:,:,numberOfFrames))



%implay(frames);