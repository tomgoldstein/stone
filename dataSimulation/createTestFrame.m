%  Create a simulated video frame.  This frame can be used for video
%  reconstruction experiments.

function frame = createTestFrame( n, time, frame )


%frame = phantom(n,n);

if ~exist('frame','var')
    frame = zeros(n,n);
end

time = min(time,n-30);

frame(end-36:end-30, 1+time:6+time) = 255;

%frame = circshift(frame,[time,time]);

end

