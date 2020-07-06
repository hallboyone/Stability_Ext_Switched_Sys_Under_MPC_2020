function [UB, gamma, sCost] = FragmentsUB(c, mode, ops)
%% SYSTEM SETUP
assumption_true = 1;

%Arrays to hold parameters
gamma = zeros(1,numel(mode)); %s.t. J*_i(x(t+1)) < gamma * J*_i(x(t))
sCost = eye(numel(mode));    %s.t. J*_i(x(t)) < S_{j,i} * J*_j(x(t))

%% BORDER CONTROL OUTER
fprintf("\nRunning the Fragments Upper Bound\n");
if(ops.run_timer)
    tic
end
%Get the set in which a switch is feasible
feas_set =  mode(1).S(end);
for m_i = 2:numel(mode)
    feas_set = feas_set & mode(m_i).S(end);
end
for m_i = 1:numel(mode)
    fun_name = ['lb',int2str(m_i)];
    feas_set.addFunction(@(x) x'*mode(m_i).P*x, fun_name);
end

for m_i = 1:numel(mode)
    %===========================================================================
    %==== Compute the largest increase in cost when switching to this mode =====
    %===========================================================================
    
    %Fragment the feas set and get the cost weights within each fragment
    [costSet, fragments, ~] = fragmentQuadCost(mode(m_i), feas_set);
    
%     figure
%     hold on
%     for i=1:numel(fragments)
%         fragments(i).addFunction(@(x) x'*costSet{i}*x, 'ub');
%         fragments(i).fplot('ub');
%     end
    for src_mode = 1:numel(mode)
        %If the src_mode is the current mode, switch cost is 1
        if src_mode == m_i
            continue;
        end
        
        
        if ops.DEBUG
            fprintf("Looking for S_cost in %d -> %d\n", src_mode, m_i);
        end
        
        %Get the terminal cost for the lower bound
        P = mode(src_mode).P;
        
        %For each fragment of the feasible region, find the largest gap
        %between the lower and upper bounds
        for frag_i = 1:numel(fragments)
            %Maximize the ratio using the Frank-Wolfe method
            maxRatio = frankWolfe(costSet{frag_i}, P, fragments(frag_i).V);
            sCost(src_mode, m_i) = max(sCost(src_mode, m_i), maxRatio);
        end
    end
    
    if ops.figs || ops.DEBUG
        %Generate plot for paper and info on the error
        filename = "boundary_error_" + m_i;
        FragmentError_2D(costSet, fragments, c.True{m_i}, c.Low{m_i},...
            ops.mesh, filename, m_i, mode(m_i).S(end));
    end
    
    %===========================================================================
    %======= Compute the smallest step down while dwelling in this mode ========
    %===========================================================================
    [costSet, fragments, Ksets] = fragmentQuadCost(mode(m_i), mode(m_i).S(end));
    for frag_i = 1:numel(fragments)
        if assumption_true
            gamma(m_i) = max(gamma(m_i),...
                1+frankWolfe(-(mode(m_i).Q+Ksets{frag_i}'*mode(m_i).R*Ksets{frag_i}), costSet{frag_i}, fragments(frag_i).V));
        else
            gamma(m_i) = max(gamma(m_i), 1+frankWolfe(-(mode(m_i).Q), costSet{frag_i}, fragments(frag_i).V));
        end
    end
end
if(ops.run_timer)
    c.time = c.time + toc;
    c.time
end

if assumption_true
    MJLSConstraints(sCost, gamma, "Fragments + Gamma Assumption");
    minDwellTimes(sCost, gamma, "Fragments + Gamma Assumption");
else
    MJLSConstraints(sCost, gamma, "Fragments");
    minDwellTimes(sCost, gamma, "Fragments");
end

UB = nan;
end