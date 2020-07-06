% THIS UPPERBOUND IS NOT CORRECT. A SIMPLE CONTER EXAMPLE CAN BE
% CONSTRUCTED. DO NOT USE!

%% SYSTEM SETUP
clearvars
clc

%Load the modes. If they are saved, they will be opened, otherwise, they
%are built.
[m, costs] = buildModes(2, 0, 1);

%Define the mesh
X = -5:0.1:5;
Y = -5:0.1:5;

c3 = repmat(struct('FreeFlow_FullPower', 0), 2, 1);

sCost = repmat(struct('FreeFlow_FullPower', ones(1,size(m, 1))), 2, 1);

gamma = repmat(struct('FreeFlow_FullPower', 0), 2, 1);

cQ = zeros(1, numel(m));
for i=1:numel(m)
    cQ(i) = min(svd(m(i).Q));
end

%% FREE ROAM - FULL POWER
d = (size(X,2) - 1)/80;
cur_d = d;
fprintf("\nRunning the Free-Roam/Full-Power Upper Bound\n");
tic
%m = buildModes(2, 1, 0);
for mode = 1:size(m, 1)

    %Get the cost of full power over the entire horizon
    U_UB = m(mode).N * maxOverSet(@(u) u'*m(mode).R*u, m(mode).U);
    [Q, ~] = BatchCost(m(mode).Q, m(mode).R, m(mode).P, m(mode).N);
    A = BatchA(m(mode).A, m(mode).N);
    Q = A'*Q*A;
    f = @(x) x'*Q*x + U_UB;
    
    %At each point in grid, compute the composite cost
    for i=1:size(X,2)
        if i > cur_d
            fprintf("=");
            cur_d = cur_d + d;
        end
        for j=1:size(Y,2)
            x = [X(i); Y(j)];
            if(m(mode).T.contains(x))
                costs{mode}.FreeFlow_FullPower(j,i) = x' * m(mode).P * x;
            else
                costs{mode}.FreeFlow_FullPower(j,i) = f(x);
            end
        end
    end
    c3(mode).FreeFlow_FullPower = max(max(costs{mode}.FreeFlow_FullPower./(costs{mode}.norm.^2)));
    gamma(mode).FreeFlow_FullPower = (1 - cQ(mode)/c3(mode).FreeFlow_FullPower);
    for k = 1:size(m, 1)
        if(k~=mode)
            sCost(mode).FreeFlow_FullPower(k) = max(max(costs{mode}.FreeFlow_FullPower./costs{k}.low));
        end
    end
end
fprintf("\n");
toc

minDwellTimes(sCost, gamma)

plotErrors(X, Y, costs)
