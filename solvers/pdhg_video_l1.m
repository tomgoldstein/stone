function [x,outs] = pdhg_video_l1( R, b, mu, order )

rows = order.n; 
cols = order.n;
numFrames = size(R,2);

%%  Allocate memory
NN = rows*cols;
x0 = zeros(NN,numFrames);
y0 = zeros(4*NN,numFrames);  %  x,y,z store the derivatives in the 3 directions


%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

%b = b*10000/max(max(abs(b))); % 10000
scale = max(max(b))-min(min(b));
b = b-mean(b(:));
b = 100*b/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma) [projectInf(y(1:end/4,:),y(end/4+1:2*end/4,:)) ; min(max(y(2*end/4+1:3*end/4,:),-1),1); min(max(y(3*end/4+1:end,:)-sigma*(R.*b),-mu),mu) ];
A = @(x) [Dx*x ; Dy*x; Dz(x) ; R.*STO(x)];
At = @(y) Dxt*y(1:end/4,:)+Dyt*y(end/4+1:2*end/4,:)+Dzt(y(2*end/4+1:3*end/4,:))+STO(R.*y(3*end/4+1:end,:));


%%  Determine timestep paramters 
opts = [];
%opts.L = .95/12;
opts.tol = .002;
%opts.tau = 0.5;  % 1.0 Works well for this problem when b has maximum size 10000
%opts.sigma = opts.L/opts.tau;
opts.maxIters = 500;  % Note: may have been 2000 for paper results
opts.verbose = true;

[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);

  
end


function dz = Dz(u)
dz = imfilter(u, [-1 1 0]);
dz(:,1) = 0;
end

function dz = Dzt(u)
dz = imfilter(u, [0 1 -1]);
dz(:,1) = -u(:,2);
dz(:,end) = u(:,end);
end










