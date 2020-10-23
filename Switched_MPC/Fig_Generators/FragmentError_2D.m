function FragmentError_2D(costSet, fragments, trueCost, lowCost, mesh, filename, mode_num, feas)
%FRAGMENTERROR_2D Generates a surf plot showing the upperbound and lower bound
%error of fragments
X_vals = mesh.value(1:mesh.numPoints);

disp("Generating error surfs...");

err_UB = nan(size(trueCost));
if(~isempty(lowCost))
    err_LB = (lowCost - trueCost)./trueCost;
    err_LB(mesh.origin) = 0;
end

%For each fragment
for f=1:numel(fragments)
    %Get the H and b that define the polyhedra Hx < b
    H = fragments(f).H(:, 1:end-1);
    b = fragments(f).H(:,end);
    
    %Get the quadratic weight of the upper bound in the current fragment
    M = costSet{f};
    
    %At each point in the mesh that has not been examined already and is in
    %the wedge
    for i = 1:mesh.numPoints
        if isnan(err_UB(i))
            x = X_vals(:,i);
            if(all(H*x<b))
                err_UB(i) = (x'*M*x - trueCost(i))/trueCost(i);
                if(err_UB(i) < -10e-8)
                    fprintf("Found smaller point in ub : %0.10f\n", err_UB(i));
                end
            end
        end
    end
end
err_UB(mesh.origin) = 0;

if mesh.dim==2
    fig = figure;
    if(isempty(lowCost)) %Just plot upperbound error
        plot(feas, 'alpha', 0)
        hold on
        surf(mesh.vals(1,:), mesh.vals(2,:), reshape(err_UB, mesh.sz));
        colorbar
        shading flat
        title("Upper Bound Error - Mode "+mode_num)
        xlabel("x_1")
        ylabel("x_2")
        fig.Position = [25 427 833 686];
    else %Plot upper and lower bound
        subplot(1,2,1)
        plot(feas, 'alpha', 0)
        hold on
        surf(mesh.vals(1,:), mesh.vals(2,:), reshape(err_UB, mesh.sz));
        colorbar
        shading flat
        title("Upper Bound Error - Mode "+mode_num)
        xlabel("x_1")
        ylabel("x_2")
        subplot(1,2,2)
        plot(feas, 'alpha', 0)
        hold on
        surf(mesh.vals(1,:), mesh.vals(2,:), reshape(err_LB, mesh.sz));
        colorbar
        shading flat
        title("Lower Bound Error - Mode "+mode_num)
        xlabel("x_1")
        ylabel("x_2")
        fig.Position = [3 408 1917 708];
    end
    %saveFig(filename, 'high', fig);
end

fprintf("Max upper bound error was %0.5f\n", max(err_UB));
if(~isempty(lowCost))
    fprintf("Max lower bound error was %0.5f\n", max(-err_LB));
end
end

