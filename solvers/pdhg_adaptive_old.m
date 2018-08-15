function [x,y,tau_vec,rp,rd,f1_vec,f2_vec ]= pdhg_adaptive(x,y,A,At,xProx, yProx,L,tol,adaptive, f1, f2,tau,sigma)

MAX_ITERATIONS = 500; %500

if ~exist('adaptive','var')
    adaptive = true;
end
if ~exist('f1','var')
    f1 = @(x,y) 0;
end
if ~exist('f2','var')
    f2 = @(x,y) 0;
end

if ~exist('tau','var')
tau = 1*sqrt(L);
end
if ~exist('sigma','var')
sigma = L/tau;
end

xr = randn(size(x));
yr = randn(size(y));
prod1 = A(xr).*yr;
prod2 = xr.*At(yr);
innerProd1 = sum(prod1(:));
innerProd2 = sum(prod2(:));

assert( ...
   abs( innerProd1 - innerProd2 )/max(abs(innerProd1),abs(innerProd2)) < 1e-8...
    ,'A and At are not formal adjoints');
    


tau_vec = zeros(MAX_ITERATIONS,1);
f1_vec = zeros(MAX_ITERATIONS,1);
f2_vec = zeros(MAX_ITERATIONS,1);
rp = zeros(MAX_ITERATIONS,1);
rd = zeros(MAX_ITERATIONS,1);

Aty = At(y);
Ax = A(x);

 %% .25/.97 or .5/.95
delta = .25;
decay = .90; %0.95, Max variation = 8.28e-6.  <delta = .5; decay = .95; p = 1; for i=1:1000; p=p*(1-delta); delta = delta*decay; end; p>

count = 0;
iter = 0;
error = 1;
while   (error > tol && iter<MAX_ITERATIONS) || iter<5
 
   iter = iter+1;
   x0 = x;
   y0 = y;
   Ax0=Ax;
   Aty0=Aty;
   tau0=tau;
   
   
   x = xProx(x-tau*Aty,tau);
   Ax = A(x);
   xh = x+(x-x0);
   Axh = 2*Ax-Ax0;
  
   
   
   y = yProx(y+sigma*Axh,sigma);
   Aty = At(y);
   
   
   dx = x-x0;
   dy = y-y0;
   
   
   r = dx/tau - (Aty-Aty0);
   d = dy/sigma - (Ax-Ax0);
     
   rp(iter) = sum(sum(abs(r)))/numel(r);
   rd(iter) = sum(sum(abs(d)))/numel(d);
   f1_vec(iter) = f1(x,y);
   f2_vec(iter) = f2(x,y);
   tau_vec(iter) = tau;
    


  if adaptive && iter>1 && max(rp(iter),rd(iter))< max(rp(iter-1),rd(iter-1))
    if rp(iter)>1.2*rd(iter)
            tau = tau/(1-delta);
            sigma = L/tau;
            delta=delta*decay;
            count = count+1;
        end
        if rp(iter)<0.8*rd(iter)
            tau = tau*(1-delta);
            sigma = L/tau;
            delta=delta*decay;
            count = count+1;
        end
  end

  primalError = rp(iter)/max(rp);
  dualError = rd(iter)/max(rd);
  error =  max(primalError,dualError);
  fprintf('%d : %f\n',iter,error);
  
end

tau_vec = tau_vec(1:iter);
f1_vec = f1_vec(1:iter);
f2_vec = f2_vec(1:iter);
rp = rp(1:iter);
rd = rd(1:iter);

count