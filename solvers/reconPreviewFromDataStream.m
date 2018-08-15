function [ preview ] = reconPreviewFromDataStream(start, window, data, ind, order )

N = order.n;
R = zeros(N*N,1);
stop = start+window-1;
R(ind(start:stop)) = 1;
b = zeros(N*N,1);
b(ind(start:stop)) = data(start:stop);
preview = makePreview(R, b, order);

end

