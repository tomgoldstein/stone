%%  Convert an image to a vector with nested dissection ordering
%   Inputs are the 2D image, and the data struct created by the method
%   "createOrderingData".

function vec = imageToNestedVector( image, order )

vec = zeros(order.n*order.n,1);
vec(order.vector) = image(:);

end

