function [x,y,taus,rp,rd ]= pdhg_01_image( R, b, mu, order )

MAX_ITERS = 500;

%% Adaptive timestep parameters
adaptive  = true;
phi = .25;
decay = .95; %0.95, Max variation = 8.28e-6.  <phi = .5; decay = .95; p = 1; for i=1:1000; p=p*(1-phi); phi = phi*decay; end; p>
Delta = 1.5;  %  must be >1

L = .9/8;
tau = sqrt(L)*3;
sigma = L/tau;




assert(size(R,2)==1,'R must be a column vector');
N = size(R,1);
assert(N==int32(N/4)*4,'N must be a power of 4');

%% Determine scale on which the recovered image lives
scale = max(b-sum(b)/sum(R)); % 0/1 transform is roughly unitary when applied random vectors with mean zero
                          %  so we'd expect mean of transformed data to be
                          %  mean of original data (i.e. the original image)


 %% Define the matrix and proximal operators
[Dx Dy] = createDifferenceOperators(order); % The derivatives we need for TV
Dxt = Dx';
Dyt = Dy';
A = @(x) [Dx*x ; Dy*x];    %  This is the graident matrix needed for PDHG
At = @(y) Dxt*y(1:end/2)+Dyt*y(end/2+1:end); %...and its transpose

rhs = mu*STO(R.*b)+mu*sum(R.*b)/sqrt(N); % Needed to form the RHS of the proximal operator equation
[xProx] = createProximals(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms
yProx = @(y,tau) min(max(y,-scale),scale);  % This re-projects y into the [-1,1] interval for PDHG


   % Intialize everything with zeros
   
taus = zeros(MAX_ITERS,1);
rp = zeros(MAX_ITERS,1);
rd = zeros(MAX_ITERS,1);
x = zeros(N,1);
y = zeros(2*N,1);
Aty = x;    
Ax = y;


count = 0;
iter = 0;
rp(1) = 1;
while  iter<1 || (max(rp(iter)/rp(1),rp(iter)/rp(1))>0.002 && iter<MAX_ITERS)
   iter = iter + 1;
   x0 = x;
   y0 = y;
   Ax0=Ax;
   Aty0=Aty;
   
   x = xProx(x-tau*Aty,tau);
   Ax = A(x);
   Axh = 2*Ax-Ax0;
   
   y = yProx(y+sigma*Axh,sigma);
   Aty = At(y);
   
   
   dx = x-x0;
   dy = y-y0;
   
   
   r = dx/tau - (Aty-Aty0);
   d = dy/sigma - (Ax-Ax0);
     
   rp(iter) = sum(sum(abs(r)))/numel(r);
   rd(iter) = sum(sum(abs(d)))/numel(d)/3;
   taus(iter) = tau;
   

  if adaptive && iter>1 && max(rp(iter),rd(iter))< max(rp(iter-1),rd(iter-1))
    if rp(iter)>rd(iter)*Delta
            tau = tau/(1-phi);
            sigma = L/tau;
            phi=phi*decay;
            count = count+1;
            [xProx] = createProximals(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms
        end
        if rp(iter)<rd(iter)/Delta
            tau = tau*(1-phi);
            sigma = L/tau;
            phi=phi*decay;
            count = count+1;
            [xProx] = createProximals(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms
        end
    end

rp(iter)

end

count

return

function xProx = createProximals(R,mu,tau,rhs)

N = size(R,1);

s = sum(R); 
sr = STO(R);
one = ones(N,1);
  % U and V are the matrices that define the rank-2 update
U = [sr one/sqrt(N)]; 
V = [mu*one/sqrt(N) mu*(sr+s/sqrt(N) )]';
  
kernel = mu*R+(1/tau);
Ainv = @(x) STO(STO(x)./kernel);    % ( mu*SRS+1/tau)^-1

invU = [Ainv(U(:,1)) Ainv(U(:,2)) ]; % compute ( mu*SRS+1/tau)^-1*U
compliment = V*invU+eye(2,2);        % Compute (I+V*Ainv*U)^-1
rank2correction = @(x) (x - invU*(compliment\(V*x))); 
Zinv = @(x) rank2correction(Ainv(x));  % ( mu*ZRZ+1/tau)^-1 

xProx = @(q,tau) Zinv(q/tau+rhs);   % min_x ||RZ x -b||^2 + (1/2tau)||x - q||^2

return
