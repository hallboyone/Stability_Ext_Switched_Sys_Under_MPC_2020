function [UB, gamma, sCost] = BoundaryPreUB(c, m, ops)
%% SYSTEM SETUP

%Arrays to hold parameters
gamma = zeros(1,numel(m)); %s.t. J^*(t + 1) < gamma * J^*(t)
c_UB = zeros(1, numel(m)); %s.t. J^UB < c_UB * ||x||^2
c_Q = zeros(1, numel(m));  %s.t. q(x,u) > c_Q * ||x||^2
sCost = ones(numel(m));    %s.t. J*_i < S_{j,i} J*_j

%Set c_Q to the smallest singular value of Q
for i=1:numel(m) 
    c_Q(i) = min(eig(m(i).Q));
end

%Cell array for the upper bound cost
UB = nan;%repmat({nan(numPoints, 1)}, numel(m), 1);

%% BORDER CONTROL FULL
fprintf("\nRunning the Pre-Set Boundary Control Upper Bound\n");
tic
for mode = 1:numel(m)
    if ops.figs || ops.DEBUG
        allFrags = [];
        allCosts = {};
    end
    for i=1:m(mode).N+1
        %Get the cell array of quad costs matrices for each fragment
        [costSet, fragments] = fragmentQuadCost(m(mode), m(mode).S(i));
        
        if ops.figs || ops.DEBUG
            allFrags = [allFrags, fragments'];
            allCosts = [allCosts, costSet];
        end
        
        %For each mode we could have switched from (i.e all but our own)
        for src = 1:numel(m)
            if src ~= mode
                %Get the P matrix for the lower bound
                P = m(src).P;
                
                %For each fragment of the feasible region, find the largest gap
                %between the lower and upper bounds
                for frag_i = 1:numel(fragments)
                    %Maximize the ratio using the Frank-Wolfe method
                    maxRatio = frankWolfe(P, costSet{frag_i}, fragments(frag_i).V);
                    sCost(src, mode) = max(sCost(src, mode), maxRatio);
                    
                    %Find the largest ratio between the cost and norm.
                    c_UB(mode) = max(c_UB(mode), frankWolfe(eye(ops.dim), costSet{frag_i}, fragments(frag_i).V));
                    gamma(mode) = max(gamma(mode), 1+frankWolfe(costSet{frag_i}, -m(mode).Q, fragments(frag_i).V));
                end
            end
        end
    end
    
    if ops.figs || ops.DEBUG
        %Figure for paper
        filename = "preset_boundary_error_" + mode;
        FragmentError_2D(allCosts, allFrags, c.True{mode}, [], ops.mesh, filename, mode, m(mode).S(end));
    end
    %Compute the gamma value
    %gamma(mode) = (1 - c_Q(mode)/c_UB(mode));
end
toc

minDwellTimes(sCost, gamma, "Pre-Set Boundary")
end