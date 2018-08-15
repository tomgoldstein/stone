% Create simulated compressive data from Shepp-Logan phantom
% Returns two 3D arrays.  Entry R(:,f) is a row selector vector for frame
% f.  Also, b(:,f) is a vector with the compressive measurements from frame
% f.  Note that this code creates multiple frames by shifting the
% shepp-logan phanton.  Each frame's data is stored in a column of R or b.
%

function [R, b] = ...
        createCompressiveMeasurements( numFrames, samplesPerFrame, order, type )

if ~exist('type','var')
    type = 'STO';
end
    
sampOrder = order.samplingOrder;
n = order.n;
N = n*n;

b = zeros(N,numFrames);
R = zeros(N,numFrames);

Z = @(x) STO(x)+sum(sum(x))/sqrt(numel(x));    %  The zero-one transform matrix

%  Loop over each frame
    for f  =1:numFrames

        % Get the indices to sample
        indicesToSample = sampOrder(1:samplesPerFrame);
        % shift the sampling position so next block starts where last one left
        % off
        sampOrder = circshift(sampOrder,-samplesPerFrame);

        % Create a synthetic test image
        image = createTestFrame(order.n,f);

        % Convert the image to a vector so we can transform it
        vec = imageToNestedVector(image,order);

        %  Apply the measurement operator
       if strcmp(type,'STO')
          vec = STO(vec);
       else
          vec = Z(vec);
       end

        % Record the row-selector matrix
        R(indicesToSample, f) = 1;

        % Record the STO modes that we measured
        b(indicesToSample, f) = vec(indicesToSample);
    end

end

