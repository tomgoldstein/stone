%  Create compressive measurements from synthetic video of moving object, 
%  and then create a low-resolution preview from the data.

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
imagesc(filtered);


