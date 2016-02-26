function [xstar, fOptim, failcount] = genericpso(fname, frestr, X0, tol, ...
                                                    maxIt, varargin)
% 	 PARTICLE SWARM OPTIMISATION
% Computes the PArticle Swarm Optimisation algorithm for fname(x),
% given the restrictions in frest(x).
% 
% [xstar, fOptim, failcount] = genericpso(fname, frestr, X0, tol, maxIt)
% 
    
if nargin < 3
    tolerance = 0.5;
    maxIter = 10;
elseif nargin < 4 
    tolerance = tol;
    maxIter = 10;
else
    tolerance = tol;
    maxIter = maxIt;
end

N = 100; % number of particles
c1 = 0.5; % velocity parameters
c2 = 1;

% Check if there are restrictions for the particles.
restrictions = ~isempty(frestr);

if restrictions == true
    % Feasibility check.
    fr = feval(frestr,X0);
    frIdx = find(fr<0);
    if ~isempty(frIdx)
        disp('Initial point not feasible.');
        xstar = [];
        fOptim = [];
        failcount=[];
        return;
    else
        % Initialize particles.
        aT = 0.5 + rand(1,N);
        x = X0*aT; % matrix
        for i=1:N % could be changed for parfor
            fr = feval(frestr,x(:,i));
            frIdx = find(fr<0);
            x(:,i) = X0;
        end
    end
else
    % Initialize particles.
    aT = 0.5 + rand(1,N);
    x = X0*aT; % matrix
end

v = zeros(size(x));

gBest = 0;
pBest = zeros(size(x));
pBestValue = zeros(N,1);
valF = zeros(N,1);

count = 0;
failc = 0;

while count <= maxIter || norm(gradBest)>tol
    parfor i=1:N
        evalRes = feval(fname, x(:,i), varargin{:});
        if pBestValue(i) < evalRes
            pBest(:,i) = x(:,i);
            pBestValue(i) = evalRes;
        end
        
    end
    [gBest, whoBest] = max(pBestValue);
    
    for i=1:N % if there are no restrictions, it might work with parfor.
        cc1 = c1*rand;
        cc2 = c2*rand;
        
        v(:,i) = v(:,i) + cc1*(pBest(:,i)-x(:,i)) + ...
                 cc2*(x(:,whoBest)-x(:,i));
        x(:,i) = x(:,i) + v(:,i);
        
        if restrictions == true
            fr = feval(frestr, x(:,i));
            frIdx = find(fr<0);
            while ~isempty(frIdx)
                failc = failc + 1;
                cc1 = cc1/2;
                cc2 = cc2/2;
                v(:,i) = v(:,i) + cc1*(pBest(:,i)-x(:,i)) + ...
                         cc2*(x(:,whoBest)-x(:,i));
                x(:,i) = x(:,i) + v(:,i);
                fr = feval(frestr, x(:,i));
                frIdx = find(fr<0);
            end
        end
    end
    count = count + 1;
    [~,gradBest] = feval(fname, x(:,whoBest), varargin{:});
    
end

xstar = x(:,whoBest);
fOptim = gBest;
failcount = failc;
                
        
        
        







































