%  Function createDifferenceOperators by Tom Goldstein
%
%  This function creates first-order difference operators
%  in the x- and y-direction.  The operators act of images that  
%  have been vectorized using the nested dissection ordering.  The 
%  returned values are sparse matrices that perform the differencing.

function [ Dx Dy ] = createDifferenceOperators( order )


rows = order.n;
cols = order.n;

%  Indexes to the location of each pixel in the vector
index = order.matrix;

% Build the derivative operator in the x (rows) direction
xInd = index; %index of each pixel
diffInd = circshift(index,[-1 0]); % index of adjacent pixels
vals = ones(rows,cols);
vals(end,:)=0;  % don't take differences along the edges
xInd = xInd(:);
diffInd=diffInd(:);
vals = vals(:);
Dx = sparse([xInd;xInd],[xInd;diffInd], [-vals;vals]); % Build sparse matrix


% Build the derivative operator in the y (rows) direction   
yInd = index;
diffInd = circshift(index,[0 -1]);
vals = ones(rows,cols);
vals(:,end)=0;  % don't take differences along the edges
yInd = yInd(:);
diffInd=diffInd(:);
vals = vals(:);
Dy = sparse([yInd;yInd],[yInd;diffInd], [-vals;vals]); % Build sparse matrix

return;
