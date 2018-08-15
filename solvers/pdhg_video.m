function [x,outs] = pdhg_video( R, b, mu, order, tvdims )

if ~exist('tvdims', 'var')
    tvdims = 2;
end

rows = order.n; 
cols = order.n;
numFrames = size(R,2);

%%  Allocate memory
NN = rows*cols;
x0 = zeros(NN,numFrames);
y0 = zeros(3*NN,numFrames);  %  x,y,z store the derivatives in the 3 directions


%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

%b = b*10000/max(max(abs(b))); % 10000
scale = max(max(b))-min(min(b));
b = 300*b/scale;
rhs = mu*STO(R.*b);

xProx =@(x,tau) STO( ...
       (1/tau+mu*R).^-1.*STO(rhs+x/tau));
yProx = @(y,sigma) [projectInf(y(1:end/3,:),y(end/3+1:2*end/3,:)) ; min(max(y(2*end/3+1:end,:),-1),1) ];
if tvdims==3
    A = @(x) [Dx*x ; Dy*x; Dz(x)];
    At = @(y) Dxt*y(1:end/3,:)+Dyt*y(end/3+1:2*end/3,:)+Dzt(y(2*end/3+1:end,:));
elseif tvdims==2
    A = @(x) [Dx*x ; Dy*x; zeros(NN,numFrames)];
    At = @(y) Dxt*y(1:end/3,:)+Dyt*y(end/3+1:2*end/3,:);
else
    assert(false,'invalid dimensionality, must be 2 or 3');
end


f1 = @(x,y) 0;
f2 = @(x,y) 0;


%%  Determine timestep paramters 
opts = [];
%opts.L = .95/12;
opts.tol = .002;
%opts.tau = 0.5;  % 1.0 Works well for this problem when b has maximum size 10000
%opts.sigma = opts.L/opts.tau;
opts.maxIters = 1500;
opts.verbose = false;

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










