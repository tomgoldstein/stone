function [ frames ] = flipud3( frames )

for f = 1:size(frames,3)
    frames(:,:,f) = flipud(frames(:,:,f));
end

end

