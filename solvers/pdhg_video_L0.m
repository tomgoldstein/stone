function [x,y,tau,rp,rd,A,At] = pdhg_video_L0( R, b, mu, order, x0,y0)

rows = order.n; 
cols = order.n;
numFrames = size(R,2);

%%  Allocate memory
NN = rows*cols;


%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

b = b*10000/max(max(abs(b))); % 10000
rhs = mu*STO(R.*b);


xProx =@(x,tau) STO( ...
       (1/tau+mu*R).^-1.*STO(rhs+x/tau));

mag =   sqrt(y0(1:end/3,:).^2+y0(end/3+1:2*end/3,:).^2)+1e-3;
   
yProx = @(y,sigma) [projectInf(y(1:end/3,:).*mag,y(end/3+1:2*end/3,:).*mag)./[mag;mag] ; min(max(y(2*end/3+1:end,:),-1),1) ];
A = @(x) [Dx*x ; Dy*x; Dz(x)];
At = @(y) Dxt*y(1:end/3,:)+Dyt*y(end/3+1:2*end/3,:)+Dzt(y(2*end/3+1:end,:));

f1 = @(x,y) 0;
f2 = @(x,y) 0;


%%  Determine timestep paramters 
L = .95/12;
tol = .01;
tau = 2.0;  % 1.0 Works well for this problem when b has maximum size 10000
sigma = L/tau;
adaptive = true;

[x,y,tau,rp,rd]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,L,tol,adaptive, f1, f2,tau,sigma);

  
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










