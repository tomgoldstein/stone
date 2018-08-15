%  Create synthetic measurements from a static scene (Shepp-Logan phantom),
%  and reconstruct from undersampled data using compressive sensing

%%  Set up the problem
 image = phantom(256,256);  % create the test image
 image(20:25,20:25) = 1;
 numberOfSamples = round(numel(image)/10);
 mu = 100;
  
%image = double(rgb2gray(imread('../images/arch.jpg')));
%numberOfSamples = round(128*128);

sigma = 0.00;
n = size(image,1);
N = n*n; % number of pixels
order = createOrderingData(n,'full'); % create the ordering to map into a vector 


   %% Create STO transform data
vec = imageToNestedVector(image, order); % vectorize the image
vec = STO(vec);       % transform it
rowsToSample = order.samplingOrder(1:numberOfSamples);

    % Record compressive measurements
b = zeros(N,1);
b(rowsToSample) = vec(rowsToSample);
b(rowsToSample) = b(rowsToSample)+randn(numberOfSamples,1)*sigma/max(b); % add noise


% create row selector matrix
R = zeros(N,1);
R(rowsToSample) = 1;
%% Call the solver to minimize L1 and recover image from sub-sampled data
[x,outs ] = pdhg_image_l1( R, b, mu, order );

rec = nestedVectorToImage(x,order);

imagesc( medfilt2(rec));