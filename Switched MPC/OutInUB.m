function [UB, gamma, sCost] = OutInUB(c, m, ops)
%% SYSTEM SETUP

%Make list of norms from 0 to the max in the state space
norms = unique(c.Norm);
numPoints = numel(norms);
%norms = linspace(0, max(c.Norm), numPoints);

%Vars to hold the turning times and costs
turningTimes = zeros(numel(m), numel(norms));
costs = zeros(numel(m), numel(norms));

%Parameters to compute and store gamma
gamma = zeros(1,numel(m));
c_UB = zeros(1, numel(m));
c_Q = zeros(1, numel(m));
for i=1:numel(m)
    c_Q(i) = min(svd(m(i).Q));
end

%Matrix for the max increase in cost when switching between modes
sCost = ones(numel(m));

UB = repmat({nan(numPoints, 1)}, numel(m), 1);

tic

%Parameters
para = struct;

%For each mode, at each norm, compute the turning time and cost ub
for mode = 1:numel(m)
    para.N = m(mode).N;
    para.lamda_Q_max = max(eig(m(mode).Q));
    para.maxSigmaA = max(svd(m(mode).A));
    para.minSigmaA = min(svd(m(mode).A));
    para.maxSigmaB = max(svd(m(mode).B));
    para.maxU = max(vecnorm(m(mode).U.V'));
    para.maxX = max(vecnorm(m(mode).S(end).V'));
    para.rho_F_max = max(vecnorm(m(mode).S(end).V'));
    para.rho_T_max = max(vecnorm(m(mode).S(1).V'));
    para.rho_T_min = min(abs(m(mode).T.H(:,end))./vecnorm(m(mode).T.H(:,1:end-1)')');
    %Terminal cost
    J_UB_T = maxOverSet(m(mode).P, m(mode).S(1));
    
    %Input cost
    J_UB_U = m(mode).N*maxOverSet(m(mode).R, m(mode).U);
    
    for i=1:numel(norms)
        alpha = norms(i);
        %[costs(mode,i), turningTimes(mode, i)] = upperBound(m(mode), norms(i));
        if alpha < para.rho_T_min
            costs(mode, i) = norm(m(mode).P, 2) * alpha^2;
            turningTimes(mode, i) = para.N;
        else
            tau = turningTimes(mode, i-1);
            while(alphaOut(alpha, tau, para) > alphaIn(para.rho_T_max, para.N-tau, para))
                tau = tau - 1;
            end
            turningTimes(mode, i) = tau;
            
            J_UB_X = 0;
            for k=0:tau
                J_UB_X = J_UB_X + alphaOut(alpha, k, para)^2;
            end
            for k=1:(m(mode).N-tau-1)
                J_UB_X = J_UB_X + alphaIn(para.rho_T_max, k, para)^2;
            end
            J_UB_X = J_UB_X * para.lamda_Q_max;
            
            costs(mode, i) = J_UB_X + J_UB_T + J_UB_U;
        end
    end
    
    
    c_UB(mode) = max(costs(mode, i)./(norms.^2));
    gamma(mode) = (1 - c_Q(mode)/c_UB(mode));
    for k = 1:numel(m)
        if(k~=mode)
            sCost(mode, k) = max(max(UB{mode}./(norms.^2*min(eig(m(k).P)))));
        end
    end
end
toc
minDwellTimes(sCost, gamma, "Out-In")
end

function m = maxOverSet(M, Set)
m = -Inf;
if(nargin < 2)
    t = 0:pi/50:2*pi;
    x = [cos(t); sin(t)];
    for i=1:size(x,1)
        m = max(fun(x(:,i)), m);
    end
else
    x = Set.V;
    for i=1:size(x,1)
        m = max(x(i,:)*M*x(i,:)', m);
    end
end
end

function [cost, tau] = upperBound(m, alpha)
%Turning time
tau = turningTime(m, alpha);

%Cost eigenvalue
lamda_Q_max = max(eig(m.Q));

%Set radius
rho_F_max = max(vecnorm(m.S(end).V'));
rho_T_max = max(vecnorm(m.S(1).V'));

% %Extreme singular values
% sigma_A_max = max(svd(m.A));
% sigma_B_max = max(svd(m.B));
% sigma_A_min = min(svd(m.A));
% 
% %Max input norm
% u_hat = max(vecnorm(m.U.V'));


%State cost
J_UB_X = 0;
for i=0:tau
    J_UB_X = J_UB_X + alphaOut(alpha, i, para)^2;
end
for i=1:(m.N-tau-1)
    J_UB_X = J_UB_X + alphaIn(rho_T_max, i, para)^2;
end
J_UB_X = J_UB_X * lamda_Q_max;

%Terminal cost
J_UB_T = maxOverSet(m.P, m.S(1));

%Input cost
J_UB_U = m.N*maxOverSet(m.R, m.U);

%Total upper bound
cost = J_UB_X + J_UB_T + J_UB_U;
end

function tau = turningTime(alpha, para)
tau = para.N;

while(alphaOut(alpha, tau, para) > alphaIn(para.rho_T_max, para.N-tau, para))
    tau = tau - 1;
end
end

function maxNorm = alphaOut(startingNorm, k, para)
%Get the max potential norm, ignoreing the feasible region constraints
maxNorm = para.maxSigmaA^k*startingNorm + para.maxSigmaB * para.maxU * sum(para.maxSigmaA.^(0:k-1));
%Apply feasible region constraints. 
maxNorm = min(para.maxX, maxNorm);
end

function maxNorm = alphaIn(endingNorm, k, para)
%Get the max potential norm, ignoreing the feasible region constraints
maxNorm = endingNorm/(para.minSigmaA^k) + (para.maxSigmaB/para.minSigmaA) * para.maxU * sum(1./para.minSigmaA.^(0:k-1));

%Apply feasible region constraints. 
maxNorm = min(para.maxX, maxNorm);
end