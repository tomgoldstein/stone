% Create a preview from undersampled data.  The method automatically
% calculates the downsampling factor by using the number of measurements.

function [preview vec]= makePreview(R, b, order )

assert(isempty(find(R~=1 & R~=0,1)),'R must contain only zero and 1');

N = size(R,1);  %  Number of pixels in high resolution
n = sum(R);     % pixels in low resolution: one measurement per pixel


delta = N/n; % down-sampling factor

%  If the measurements are complete, then we just take inverse transform
if delta==1
    vec = STO(b);
    preview = nestedVectorToImage(vec, order);
    return;
end

%  Average blocks of consecutive measurements.  Blocks have length delta
blockData = reshape(b,[delta, N/delta]);
blockR = reshape(R,[delta, N/delta]);
means = sum(blockR.*blockData)./sum(blockR);

% Take the averages, and up-sample them back to the original resolution.
% Now, each length delta block of vec has constant value, and this value is
% equal to the mean of the samples in the corresponding block of 'b'.
blocks = kron(means,ones(delta,1));
vec = blocks(:);

vec = STO(vec);

preview = nestedVectorToImage(vec, order);

%  Display the preview at the low resolution by down-sampling.  This step
%  is optional.  Without it, you'd get the low res preview, but displayed
%  at a higher resolution.
rootDelta = round(sqrt(delta));
preview = preview(1:rootDelta:end,1:rootDelta:end);

end

