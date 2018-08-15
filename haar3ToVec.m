% Map a wavelet transform object into a vector

function [ rval ] = haar3ToVec( H )
dec = H.dec
rval = []
for i = 1:numel(dec)
    a = dec{i};
    rval = [rval ; a(:)];
end

end

