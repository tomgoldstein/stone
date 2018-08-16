%  Create compressive measurements from synthetic video of moving object, 
%  and then create a low-resolution preview from the data.

%% Make preview from Shepp-Logan phantom
n = 256; % use n by n image
N = n*n; % number of pixels
subsample_factor = 4; % How much to sub-sample data for the preview.  subsample_factor=4 means use 1/4 the data.
order = createOrderingData(n,'semi'); % create the ordering to map into a vector 
image = phantom(n);

vector = imageToNestedVector(image,order);
b = STO(vector);  %  Obtain measurements
R = zeros(N,1);    %  a mask to record which data we've measured

% randomly sub-sample data
R(1:subsample_factor:N) = 1;  % measure a quarter of the data, and mark the locations of measured data in the mask

% Create preview
[preview vec] = makePreview(R,b,order);
filtered = medfilt2(preview);
figure;
subplot(1,2,1);
imagesc(filtered);
title('Shepp-Logan (25%%) data');


%%  Set up the problem
n = 256; % use n by n image
N = n*n; % number of pixels
order = createOrderingData(n,'semi'); % create the ordering to map into a vector 

np = 128; % preview resolution


%%  Create synthetic compressive measurements for the frame
%   'Frames' is a cell array. Each entry is a struct with a row-selector
%   and vector of measurements for a frame
numberOfFrames = 4;
samplesPerFrame = np*np/numberOfFrames;
[R b] = createCompressiveMeasurements(numberOfFrames, samplesPerFrame, order, 'STO');

% Average all frame data together
R = sum(R,2);  
b = sum(b,2);

%% Create preview
[preview vec] = makePreview(R,b,order);
filtered = medfilt2(preview);
subplot(1,2,2);
imagesc(filtered);
title('Synthetic moving object (25%% data)');


