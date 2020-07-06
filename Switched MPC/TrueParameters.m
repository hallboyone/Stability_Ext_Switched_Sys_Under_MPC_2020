function [gamma, sCost] = TrueParameters(c, m, ops)
%Arrays to hold parameters
gamma = ones(1,numel(m)); %s.t. J^*(t + 1) < gamma * J^*(t)
sCost = ones(numel(m));    %s.t. J*_i < S_{j,i} J*_j

numPoints = ops.mesh.numPoints;

%Make matrix of all X_vals. Greatly increases speed
X_vals = ops.mesh.value(1:numPoints);

for mode=1:numel(m)
    %Get the minimum value of gamma
    Q = m(mode).Q;
    R = m(mode).R;
    A = m(mode).A;
    B = m(mode).B;
    for i = 1:numel(c.Valid{mode})
        idx = c.Valid{mode}(i);
        x = X_vals(:,idx);
        u = c.Inputs{mode}(:,idx);
%         u = cell2mat( m(mode).MPC_Cntrl(x));
%         u = u(:,1);
        gamma(mode) = min(gamma(mode), (x'*Q*x + u'*R*u)/(c.True{mode}(idx)), 'omitnan');
        if ops.DEBUG
            if(norm(x,2) > 10e-2)
                if(1-gamma(mode) < m(mode).MPCCost(A*x+B*u)/c.True{mode}(idx))
                    disp(-1+gamma(mode) + m(mode).MPCCost(A*x+B*u)/c.True{mode}(idx))
                end
            end
        end
    end
    gamma(mode) = 1 - gamma(mode);
    
    
    for src = 1:numel(m)
        if (src ~= mode)
            sCost(src, mode) = max(c.True{mode}./c.True{src});
        end
    end
end

minDwellTimes(sCost, gamma, "True Cost")
end