% Convert a nested-dissection ordered vector into a 2D image

function image = nestedVectorToImage( vec, order )

columnMajor = vec(order.vector);

image = reshape(columnMajor, [order.n, order.n]);


end

