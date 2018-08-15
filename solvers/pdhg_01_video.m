%% NOTE:  NEED TO FIX THE METHOD TO BUILD PROXIMALS
function [u,rp ] = pdhg_01_video( R, b, mu, order )

MAX_ITERS = 500;
numFrames = size(R,2);

%% Adaptive timestep parameters
adaptive  = true;
phi = .25;
decay = .95; %0.95, Max variation = 8.28e-6.  <phi = .5; decay = .95; p = 1; for i=1:1000; p=p*(1-phi); phi = phi*decay; end; p>
Delta = 3;  %  must be >1

L = .9/12;
tau = .1;
sigma = L/tau;


%%  Allocate memory
numFrames = size(R,2);
N = size(R,1);
assert(N==int32(N/4)*4,'N must be a power of 4');

%% Determine scale on which the recovered image lives
scale = mean(max( bsxfun(@minus,b,(sum(b)./sum(R)))  )) ; % 0/1 transform is roughly unitary when applied random vectors with mean zero
                          %  so we'd expect mean of transformed data to be
                          %  mean of original data (i.e. the original image)
spaceThreshold = scale;
timeThreshold = scale;
                          
%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

rhs = bsxfun(@plus,mu*STO(R.*b),mu*sum(R.*b)/sqrt(N) ); % Needed to form the RHS of the proximal operator equation

xProx = createMultiFrameProximal(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms

%% Allocate Memory
u=zeros(N,numFrames);
x=zeros(N,numFrames);  %  x,y,z store the derivatives in the 3 directions
y=zeros(N,numFrames);
z=zeros(N,numFrames);

rp = 1;
count=0;
iter = 1;
while max(rp(end)/rp(1))>0.002 && iter<500
%while iter<50 

    u0 = u;
    x0 = x;
    y0 = y;
    z0 = z;

    ubar = u-tau*(Dxt*x+Dyt*y+Dzt(z));
    
    for f=1:numFrames
        u(:,f) = xProx{f}(ubar(:,f),tau);
    end  
    uh = u+(u-u0);
  
    
    x = x+sigma*(Dx*uh);
    y = y+sigma*(Dy*uh);
    z = z+sigma*(Dz(uh));
    normalize = sqrt(x.*x+y.*y);
    normalize = max(normalize,spaceThreshold);
    x = x./normalize*spaceThreshold;
    y = y./normalize*spaceThreshold;
    z = max(min(z,timeThreshold),-timeThreshold);
   
    
  
   
   du = (u-u0)/tau;
   dx = (x-x0)/sigma;
   dy = (y-y0)/sigma;
   dz = (z-z0)/sigma;  
   
   rp(iter) = sum(sum(abs(du)))/numel(du);
   rd(iter) = sum(sum(abs(dx)))/numel(dx)+sum(sum(abs(dy)))/numel(dy)+sum(sum(abs(dz)))/numel(dz);
   taus(iter) = tau;
   

  if adaptive && iter>1 && max(rp(iter),rd(iter))< max(rp(iter-1),rd(iter-1))
    if rp(iter)>rd(iter)*Delta
            tau = tau/(1-phi);
            sigma = L/tau;
            phi=phi*decay;
            count = count+1;
           xProx = createMultiFrameProximal(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms
        end
        if rp(iter)<rd(iter)/Delta
            tau = tau*(1-phi);
            sigma = L/tau;
            phi=phi*decay;
            count = count+1;
            xProx = createMultiFrameProximal(R,mu,tau,rhs); %  The proximal operators for the l1 and l2 terms
        end
    end

   tau
    disp([num2str(iter) ' : ' num2str(rp(end)/rp(1))]);
    
    iter = iter+1;
    
end
 count
  
return


function dz = Dz(u)
dz = imfilter(u, [-1 1 0]);
dz(:,1) = 0;
return

function dz = Dzt(u)
dz = imfilter(u, [0 1 -1]);
dz(:,1) = -u(:,2);
dz(:,end) = u(:,end);
return



function xProx = createMultiFrameProximal(R,mu,tau,rhs)

nFrames = size(R,2);
xProx = {};
for f=1:nFrames
    xProx{f} = createSingleProximal(R(:,f),mu,tau,rhs(:,f));
end

return

function xProx = createSingleProximal(R,mu,tau,rhs)

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






