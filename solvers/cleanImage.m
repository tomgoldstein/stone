function [ frames ] = cleanImage( frames, filter )
small = min(frames(:));
big = max(frames(:))*1.1;
frames = (frames-small)*255/(big-small);


if ~exist('filter','var') || filter
    av = mean(frames(:));
    frames = frames-av;

    for i=1:size(frames,3)
        frames(:,:,i) = medfilt2(frames(:,:,i));
    end
    frames = medfilt3(frames); 

    frames = frames+av;
end

end
