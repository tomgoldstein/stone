%  Open a directory of png images, and create a stream of measurements from
%  them.

function [ data, finalImage ] = createMeasurementsFromImages( directory, numFrames, samplesPerFrame, order, addAnomoly, measurementCycleLength)

if ~exist('addAnomoly','var')
    addAnomoly  = false;
end

if ~exist('measurementCycleLength','var')
    measurementCycleLength  = length(order.samplingOrder);
end

data = zeros(numFrames*samplesPerFrame,1);
files = dir([directory '/*.png']);

count = 1;
for f = 1:numFrames
    im = rgb2gray(imread([directory '/'  files(f).name]));
    im = double(im);
    if f>numFrames/2 && addAnomoly
        im(end/2:end/2+4,end/2:end/2+4) = 2e3; % A bright anomaly
    end
    
    vec = imageToNestedVector(im,order);
    vec = STO(vec);
    
    for s=1:samplesPerFrame
        index = mod(count-1,measurementCycleLength)+1;
        index = order.samplingOrder(index);
        data(count) = vec(index);
        count = count+1;
    end

end

finalImage = im;


end

