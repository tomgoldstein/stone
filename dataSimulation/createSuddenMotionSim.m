%  Open a directory of png images, and create a stream of measurements from
%  them.

function [ data, ind ] = createSuddenMotionSim(dataBeforeMove, dataAfterMove, dataPerPixelMove, order)

numData = dataBeforeMove+dataAfterMove;
n = order.n;
data = zeros(numData,1);
ind = zeros(numData,1);

cam = double(imread('../images/cameraman.tif'));
%cam = cam(1:2:end,1:2:end);

shift = 0;
im = createTestFrame(n,0,cam);   
vec = imageToNestedVector(im,order);
vec = STO(vec);
for count = 1:numData
    if count>=dataBeforeMove && mod(count-dataBeforeMove,dataPerPixelMove)==0;
        shift = shift+1;
        im = createTestFrame(n,shift, cam);
        vec = imageToNestedVector(im,order);
        vec = STO(vec);
       % disp('New Frame')
    end
    index = mod(count-1,n*n)+1;
    index = order.samplingOrder(index);
    data(count) = vec(index);
    ind(count) = index;
end

end

