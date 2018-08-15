%%      pdhg_adaptive.m  by Tom Goldstein
%   This method solves the saddle-point problem
%
%    max_y min_x f(x) + <Ax,y> - g(y)
%
%   Using an adaptive PDHG method.  The required inputs are
%   "x" - An Nx1 column vector with the initial guess for x
%   "y" - An Mx1 column vector with the initial guess for y
%   "A" - Function handle to an MxN linear operator
%   "At"- Function handle to the transpose (Hermitian) of A
%   "xProx"- A function handle of the form @(x,tau)<do stuff>.
%       This function computes:  min_x  f(x)+(1/2/tau)||x-x0||^2 
%   "yProx"- A function handle of the form @(y,sigma)<do stuff>.
%       This function computes:  min_y  g(y)+(1/2/sigma)||y-y0||^2 
%   "opts" - an optional struct containing various options that the user 
%       may choose to set.  The options are described below in the function
%       "setDefaults" at the bottom of this file
%   Note: This code can handle complex matrices provided A and At are
%   Hermitian adjoints.

function [x,outs]= pdhg_adaptive(x,y,A,At,xProx, yProx,opts)

%% Make sure 'opts' struct is filled with all the options
if ~exist('opts','var')
    opts = [];
end
opts = setDefaults(opts,x,A,At);
testAdjoints(A,At,x,y);


%% Get some commonly used values from the 'opts' struct
tau = opts.tau;
sigma = opts.sigma;
maxIters = opts.maxIters;
a = opts.a;
L = opts.L;
Delta = opts.Delta;

%% Allocate space for the returned variables in the 'outs' struct
outs = [];
outs.tau = zeros(maxIters,1);
outs.sigma = zeros(maxIters,1);
outs.f1 = zeros(maxIters,1);
outs.f2 = zeros(maxIters,1);
outs.p = zeros(maxIters,1);
outs.d = zeros(maxIters,1);

%% Initialize some values
updates = 0;
Ax = A(x);
Aty = At(y);

%% Begin Iteration
for iter = 1:maxIters
   
   % store old iterates
   x0 = x;
   y0 = y;
   Ax0=Ax;
   Aty0=Aty;
    
   % primal update
   x = xProx(x-tau*Aty,tau);
   Ax = A(x);
   Axh = 2*Ax-Ax0;
    
   % dual update
   y = yProx(y+sigma*Axh,sigma);
   Aty = At(y);
     
   % compute and store residuals
   dx = x-x0;
   dy = y-y0;
   r = dx/tau - (Aty-Aty0);
   d = dy/sigma - (Ax-Ax0);
   primal = norm(r(:),1)/numel(r);
   dual = norm(d(:),1)/numel(d);
   outs.p(iter) = primal;
   outs.d(iter) = dual;
   
   if opts.verbose
    fprintf('%d: residual = %d,   change = %d\n',iter,max(primal,dual)/max(outs.p(1),outs.d(1)), norm(x-x0,'fro') );
   end
   
   % store various values that we wish to track
   outs.f1(iter) = opts.f1(x,y,x0,y0,tau,sigma);
   outs.f2(iter) = opts.f2(x,y,x0,y0,tau,sigma);
   outs.tau(iter) = tau;
   outs.sigma(iter) = sigma;
   
   % Check stopping condition
   if isfield(opts,'eps')
       stopNow = max(primal,dual)<opts.eps || iter == opts.maxIters;
   else
       stopNow = max(primal,dual)/max(outs.p(1),outs.d(1)) < opts.tol...
                                            || iter == opts.maxIters;
   end   
   if iter>5 && stopNow
     outs.y = y;
     outs.p = outs.p(1:iter);
     outs.d = outs.d(1:iter);
     outs.f1 = outs.f1(1:iter);
     outs.f2 = outs.f2(1:iter);
     outs.updates  = updates;
     outs.tau = outs.tau(1:iter);
     outs.sigma = outs.sigma(1:iter);
     outs.iters = iter;
    return;
   end
     
  % Test the backtracking/stability condition
  Axy = 2*real(sum(sum((Ax-Ax0).*conj(dy))));
  Hnorm = norm(dx,'fro')^2/tau+ norm(dy,'fro')^2/sigma;
  if opts.backtrack && opts.gamma*Hnorm<Axy
      x=x0;
      y=y0;
      Ax=Ax0;
      Aty=Aty0;  
      decay = opts.b*opts.gamma*Hnorm/Axy;
      tau = tau*decay;
      sigma = sigma*decay;
      L = L*decay*decay; 
  end
  
  %  Perform adaptive update
  if opts.adaptive && iter>1 && max(primal,dual)< max(outs.p(iter-1),outs.d(iter-1))  
    if  primal>Delta*dual
            tau = tau/(1-a);
            sigma = L/tau;
            a=a*opts.eta;
            updates = updates+1;
    end
    if  primal < dual/Delta
            tau = tau*(1-a);
            sigma = L/tau;
            a=a*opts.eta;
            updates = updates+1;
    end
  end
  
   
end

function testAdjoints(A, At, x, y)
%  Test that A and At are adjoint of one another
xr = randn(size(x));
yr = randn(size(y));
prod1 = A(xr).*yr;
prod2 = xr.*At(yr);
innerProd1 = sum(prod1(:));
innerProd2 = sum(prod2(:));

assert( ...
   abs( innerProd1 - innerProd2 )/max(abs(innerProd1),abs(innerProd2)) < 1e-8...
    ,'A and At are not formal adjoints');
    
return


%% Fill in the struct of options with the default values
function opts = setDefaults(opts,x0,A,At)



%  verbose:  Determines whether converge data is printed on each step
if ~isfield(opts,'verbose') 
    opts.verbose = false;
end

%  L:  The reciprocal of the spectral radius of A'A.
%  Approximate the spectral radius of A'A if we don't know L
if ~isfield(opts,'L') || opts.L<=0
    x = randn(size(x0));
    transform = At(A(x));
    specRadius = norm(transform,'fro')/norm(x,'fro');
    opts.L = 2/specRadius;
end

%  maxIters: The maximum number of iterations
if ~isfield(opts,'maxIters')
    opts.maxIters = 1000;
end
% tol:  The relative decrease in the residuals before the method stops
if ~isfield(opts,'tol') % Stopping tolerance
    opts.tol = 1e-2;
end
% adaptive:  If 'true' then use adaptive method.
if ~isfield(opts,'adaptive')    %  is Adaptive?
    opts.adaptive = true;
end

% backtrack:  If 'true' then use backtracking method.
if ~isfield(opts,'backtrack')    %  is backtracking?
    opts.backtrack = true;
end

% f1:  An optional function that is computed and stored after every
% iteration
if ~isfield(opts,'f1')          % These functions get evauated on each iterations, and results are stored
    opts.f1 = @(x,y,x0,y0,tau,sigma) 0;
end
% f2:  An optional function that is computed and stored after every
% iteration
if ~isfield(opts,'f2')
    opts.f2 = @(x,y,x0,y0,tau,sigma) 0;
end
% tau:  The intial stepsize for the primal variables
if ~isfield(opts,'tau')         % starting value of tau
    opts.tau = sqrt(opts.L);
end
% sigma: The intial stepsize for the dual variables
if ~isfield(opts,'sigma')       % starting value of sigma
    opts.sigma = opts.L/opts.tau;
end

%% Adaptivity parameters
if ~isfield(opts,'a')   %  Intial adaptive update strength for stepsizes
    opts.a = .5;
end
if ~isfield(opts,'eta') %  How fast does the adaptivity level decay
    opts.eta = .95;
end
if ~isfield(opts,'Delta') % update stepsizes when primal/dual ratio exceeds Delta
    opts.Delta = 1.5;
end
if ~isfield(opts,'gamma') % Used to determine when need to backtrack to maintain positivity conditions
    opts.gamma = .75;
end
if ~isfield(opts,'b')  % Adaptivity parameter used for backtracking update
    opts.b = .95;
end

return

