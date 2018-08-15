%%    createOrderingData(n) by Tom Goldstein
%  This method helps you map a 2D array of pixels into a 1D vector
%  using the nested dissection ordering.  It returns a 2D array of
%  integers, and an integer vector.  The number at each location of 
%  the 2D array is the vector index that the element gets mapped to.
%    The method also produces 'samplingOrder', which is a vector of 
%  row indices.  This vector contains the order in which the STO
%  transform rows must be sample in order to allow previews at all
%  resolutions.
%    The vector is ordering that maps a column-major vector into a 
%  nested dissection order.
%    The input argument is the side length of the nxn image.
%
% "type" is a string, either "full" or "semi"

function order = createOrderingData(n, type)


filename = ['precomputed/' num2str(n) '_' type '.mat' ];

location = mfilename('fullpath');
location = location(1:end-18);
filename = [location filename];

if exist(filename, 'file')==2
    order = load(filename);
    order = order.order;
    return;
end


assert(strcmp(type,'full') ||strcmp(type,'semi'), 'Ordering type must be either full, for full random, or semi, for semi random');

% Seed the random number generator
%rand('seed',0);
rng('default');
rng(0);
a = zeros(n,n);
a(1,1) = 1;
a(n,n) = n*n;

assert(n == 2^int32(log(n)/log(2)),'Size must be a power of 2');

%%  Create 2D array with nested dissection indices
gap = n;
while gap>1 

    for r = 1:gap:n;
        for c=1:gap:n;
            low = a(r,c);
            high = a(r+gap-1,c+gap-1);
            inc = (high-low+1)/4;
            
            lows = [low, low+inc, low+2*inc, low+3*inc];
            highs = [low+inc-1,low+2*inc-1,low+3*inc-1,low+4*inc-1];
            
            perm = randperm(4);
            %perm = [1 2 3 4];
            lows = lows(perm);
            highs = highs(perm);
            
            
            a(r,c) = lows(1);
            a(r+gap/2-1,c+gap/2-1) = highs(1);
           
            a(r,c+gap/2) = lows(2);
            a(r+gap/2-1,c+gap-1) = highs(2);
            
            a(r+gap/2,c+gap/2) = lows(3);
            a(r+gap-1,c+gap-1) = highs(3);
            
            a(r+gap/2,c) = lows(4);
            a(r+gap-1,c+gap/2-1) = highs(4);
                  
            
        end
    end
    
    gap = gap/2;
end


%% Create the ordering of row indices for sampling
ind = (1:n*n)';
blockSize = 1;
while blockSize<n*n
    spot = 1;
    while spot<n*n
       chunk = ind(spot:spot+4*blockSize-1);
       ordering = [1:blockSize ;...
           blockSize+1:2*blockSize;...
           2*blockSize+1:3*blockSize;...
           3*blockSize+1:4*blockSize];
       
       perm = 1:4;
       if strcmp(type,'full') 
           perm = randperm(4);
       end
       
%        seed = mod(spot*16/n/n,16);
%        %seed = spot;
%        %rng(seed);
%        perm = randperm(4);
       
       
       ordering(perm,:) = ordering;
       ordering = ordering(:);
      
       chunk = chunk(ordering); 
       
       ind(spot:spot+4*blockSize-1) = chunk;
       
       spot = spot+4*blockSize;
    end
    
    blockSize = blockSize*4;
end



order = [];
order.n = n;
order.matrix = a;
order.vector = a(:);
order.samplingOrder = ind;


save(filename,'order');


return;

