function [x,outs ] = pdhg_image_l1( R, b, mu, order, guess )

rows = order.n; 
cols = order.n;


%%  Allocate memory
NN = rows*cols;
x0=zeros(NN,1);
y0=zeros(3*NN,1);


%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

b = b*100/max(max(abs(b)));

xProx =@(x,tau) x;
yProx = @(y,sigma)[projectInf(y(1:end/3,:),y(end/3+1:2*end/3,:) ) ; min(max(y(2*end/3+1:end,:)-sigma*(R.*b),-mu),mu) ];
A = @(x) [Dx*x ; Dy*x; R.*STO(x)];
At = @(y) Dxt*y(1:end/3,:)+Dyt*y(end/3+1:2*end/3,:)+STO(R.*y(2*end/3+1:end,:)) ;

%  TVL1
% xProx =@(x,tau) x;
% yProx = @(y,sigma)[projectInf(y(1:end/3,:),y(end/3+1:2*end/3,:) ) ; min(max(y(2*end/3+1:end,:)-sigma*b,-mu),mu) ];
% A = @(x) [Dx*x ; Dy*x; x];
% At = @(y) Dxt*y(1:end/3,:)+Dyt*y(end/3+1:2*end/3,:)+y(2*end/3+1:end,:) ;


%%  Determine timestep paramters 
opts = [];
%opts.L = 0.95/8;
opts.tol = .01;
opts.maxIters = 500;
opts.verbose = true; 

if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0); sign(R.*STO(x0))];
end

[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);

  


% %en = 0;
% rp = 1;
% rd = 1;
% 
% 
% rhs = mu*STO(R.*b);
% 
% iter = 1;
% while max(rp(end)/rp(1),rd(end)/rd(1))>0.01 && iter<500
% %while iter<50 
% 
%     u0 = u;
%     x0 = x;
%     y0 = y;
%     
%     uhat = u-tau*(Dxt*x+Dyt*y);
%     
%     u = rhs+uhat/tau;
%     u = STO(u);
%     u = u./(mu*R+1/tau);
%     u = STO(u);
%     
%   
%     uh = u+(u-u0);
%   
%     
%     x = x+sigma*Dx*uh;
%     y = y+sigma*Dy*uh;
%     normalize = sqrt(x.*x+y.*y);
%     normalize = max(normalize,1);
%     x = x./normalize;
%     y = y./normalize;
%     
%    
%     rp(iter) = norm(u-u0,'fro')/tau;
%     rd(iter) = sqrt(norm(x-x0,'fro')^2+norm(y-y0,'fro')^2)/sigma/2.0;
%     %en(iter) = sum(sum( sqrt((Dx*u).^2+(Dy*u).^2)))+mu/2*norm(R.*hadamardTransform(B.*u)/scale-s,'fro')^2;
%    
%    
%     disp([num2str(iter) ' : ' num2str(rp(end))]);
%     
%     iter = iter+1;
%     
%     
% end
% 
% image = nestedVectorToImage(u, data);
%     
% return
