%  Reconstruct DMD data as a single frame (i.e., not a multiframe video).  
%  First, use a direct reconstruction by applying the transform to all of 
%  the data.  If the data does not represent a static scene the result will 
%  have a lot of motion aliasing and look noisy.  If the true scene does 
%  not move, then the result should look good.
%     Next, create a preview from the data.


load 'dmd_measurements/car/256/32x/data1.mat';

Nr = 256; % DMD resolution
%Nr = sqrt(length(data));
Nprev = 32; % preview resultion
Mcs = round(Nr*Nr/10); % preview resultion

Nd = size(data,1); % number of data
Np = Nr*Nr;  % number of pixels

%%  re-order the data
order = createOrderingData(Nr,'full');
b  = zeros(Np,1);
b0(order.samplingOrder) = data(1:Nr*Nr);

%% Do a conventional forward reconstruction
vec = STO(b0);

figure;
% Display result
subplot(2,2,1);
im = nestedVectorToImage(vec,order);
imagesc(im);
title('Full Resolution');
drawnow;

%% Do a preview

% sub-sample the data
b  = zeros(Np,1);
R  = zeros(Np,1);

rowsToSample = order.samplingOrder(1:Nprev*Nprev);


b(rowsToSample) = data(1:Nprev*Nprev);
R(rowsToSample) = 1;

%  get the preview
[preview vec] = makePreview(R,b,order);

% Display result
subplot(2,2,2);
im = nestedVectorToImage(vec,order);
im = cleanImage(im);
imagesc(im);
title('preview');
drawnow;

%% Sub-sampled compressive
rowsToSample = order.samplingOrder(1:Mcs);
b(rowsToSample) = data(1:Mcs);
R(rowsToSample) = 1;
subplot(2,2,3);
mu = 10;
[x,outs ] = pdhg_image( R, b, mu, order );
rec = nestedVectorToImage(x,order);
rec = cleanImage(rec);
imagesc(rec);
title('Compressive_l2');
drawnow;

subplot(2,2,4);
mu = 5;
[x,outs ] = pdhg_image_l1( R, b, mu, order );
rec = nestedVectorToImage(x,order);
rec = cleanImage(rec);
imagesc(rec);
title('Compressive_l1');
drawnow;

colormap gray;