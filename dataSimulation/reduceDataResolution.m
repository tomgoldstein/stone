%  args:  a set of frame data, one frame per column, or compressive data
%  output:  compressive measurements of each frame at a lower resolution

function [ r,b, order ] = reduceStreamResolution( R, B ,Order )

numFrames = size(R,2);
numPix = size(R,1);

r = zeros(numPix/4,numFrames);
b = zeros(numPix/4,numFrames);


for frameNumber = 1:numFrames
   Rf = R(:,frameNumber);
   Bf = B(:,frameNumber);
   
   rf = reshape(Rf,[4,numPix/4]);
   bf = reshape(Bf,[4,numPix/4]);
   
   rfs = sum(rf,1);
   bfs = sum(bf,1);
   
   
   bf = bfs./(rfs+(rfs==0));
   rf = rfs;

   b(:,frameNumber) = bf;
   r(:,frameNumber) = rf;
  
end


order = [];
order.n = Order.n/2;
order.matrix = ceil(Order.matrix(1:2:end,1:2:end)/4);
order.vector = order.matrix(:);
order.samplingOrder = ceil(Order.samplingOrder/4);
