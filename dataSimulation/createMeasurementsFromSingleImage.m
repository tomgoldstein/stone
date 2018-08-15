%  Takes in a image and a number of measurements.  This method takes the
%  first Nm measurements in the sampling order stored in 'order', and then
%  outputs a row selector, and a data vector.

function [ R,b ] = createMeasurementsFromSingleImage( im, Nm, order )


Np = order.n^2; %  number of pixels
R = zeros(Np,1);
b = zeros(Np,1);

samples = order.samplingOrder(1:Nm);

R(samples) = 1;

vec = imageToNestedVector(im,order);
vec = STO(vec);

b(samples) = vec(samples);

end

