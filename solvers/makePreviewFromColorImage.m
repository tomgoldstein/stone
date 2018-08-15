%  Create under-sampled random STO measurements from the channels of a
%  color image.  Then create a preview.  The number of measurements is 
%  (number of pixels)/downFactor^2.  Gaussian noise of standard deviation
%  sigma is added to the measurements before reconstruction.
%     NOTE:  Input image must be scaled from 0-1

function [ preview ] = makePreviewFromColorImage( im, downFactor, sigma, order )


numChannels = size(im,3);

%%  Set up the problem
n = size(im,1); % use n by n image
N = n*n; % number of pixels


%%  Create synthetic compressive measurements for the frame
%   'Frames' is a cell array. Each entry is a struct with a row-selector
%   and vector of measurements for a frame

np = n/downFactor;
numSamples = np*np;
R = zeros(N,1);
R(order.samplingOrder(1:numSamples))=1;

preview = im;

for c = 1:numChannels

b = imageToNestedVector(im(:,:,c),order);
b = STO(b)+sigma*randn(N,1);

b = R.*b;

%  Get separate previews for each channel
p = makePreview(R,b,order);
preview(:,:,c) = p;

end

preview = max(preview,0);
preview = min(preview,1);

end

