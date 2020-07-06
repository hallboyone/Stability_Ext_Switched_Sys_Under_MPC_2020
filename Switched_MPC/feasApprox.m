function Xf = feasApprox(sys, vects)
dim = size(sys.A, 1);

Xf = Polyhedron([eye(dim); -eye(dim)], .01 * ones(2*dim, 1));

active_V = Xf.V;


placed_V = [];

for i=1:vects
    %make room for the new set [x1, x2, ...] ==> [x1, 0, x2, 0,...
    
    %For each active vertex, conduct a binary search till within tol
    for k=1:size(active_V, 1)
        delta = 1;
        v = active_V(k, :); %(placed_V(2*k-1,:) - placed_V(mod(2*k,placed_V)+1,:))/2;
        while(delta >= 1/256)
            %Scale the active point
            v = v * (1 + delta);
            %Check if the scaled vertex is valid
            if(any(isnan(cell2mat(sys.MPC_Cntrl(v'))))) %Infeasible case
                v = v / (1 + delta);
                delta = delta/2;
            end
            
        end
        active_V(k, :) = v;
    end
    if(i==1)
        placed_V = active_V;
    else
        old_placed_V = placed_V;
        placed_V = zeros(size(placed_V, 1)*2, dim);
        placed_V(2*(1:size(old_placed_V, 1))-1, :) = old_placed_V;
        placed_V(placed_V(:,1)==0, :) = active_V;
    end
    active_V = (placed_V + placed_V([2:end, 1], :))/2;
end

Xf = Polyhedron(placed_V);
end

function avg_pts = getAvgOfFacets(P)
P.minHRep

end