m = buildModes(2, 0, 1);
m = m(1);

T = m.T;

tic
pts = stepOver(T);
toc

plot(T)
hold on
plot(pts(1, :), pts(2, :), '*')
x = sdpvar(2,1);
u = sdpvar(1,1);
objective = u'*m.R*u + (m.A*x+m.B*u)'*m.P*(m.A*x+m.B*u);
constraints = [m.T.H(:, 1:end-1)*(m.A*x+m.B*u) <= m.T.H(:, end),...
               m.U.H(:, 1:end-1)*u <= m.U.H(:, end)];

ops = sdpsettings('verbose', 0, 'solver', 'mosek');
c = optimizer(constraints, objective, ops, x, u);

K = [c(pts(:,end)), m.K * T.V(7,:)', m.K * T.V(3,:)']*inv([[pts(:,end);1], [T.V(7,:)';1], [T.V(3,:)';1]]);
s = K(end);
K = K(1:end-1);
X_LQR1 = Polyhedron([K; -K],[1-s; 1+s]) & m.X;
plot(X_LQR1)
hold on
plot(m.T)
plot(Polyhedron(((m.A + m.B*K)/(T.V' - m.B*s))'))

function pts = stepOver(PolySet)
eps = 0.05;
dim = PolySet(1).Dim;
lns_per_facet = ((dim-1)^2 + (dim-1))/2;
ln_idx = zeros(2, lns_per_facet);

ln_num = 1;
for i=2:dim
    for j=1:i-1
        ln_idx(:, ln_num) = [j;i];
        ln_num = ln_num + 1;
    end
end

num_edges = 0;
for i=1:numel(PolySet)
    PolySet(i).minVRep();
    PolySet(i).minHRep();
    if(PolySet(i).Dim ~= dim)
        error("Polyhedron have different diminsions!");
    end
    num_edges = num_edges + size(PolySet(i).H, 1);
end

pts = nan * zeros(dim, num_edges);
cur_pt_idx = 1;
M = zeros(lns_per_facet, dim);
for i=1:numel(PolySet)
    [I, ~] = find((PolySet(i).H(:,1:end-1) * 1.000001*PolySet(i).V' > PolySet(i).H(:, end))');
    I = reshape(I, [2,size(I,1)/2]);
    
    %For each facet of the polyhedron, get a point just beyond it
    for idx=1:size(I,2)
       v =  PolySet(i).V(I(:,idx), :);
       avg_v = sum(v)/dim;
       for k=1:lns_per_facet
           M(k, :) = [1, -1] * v(ln_idx(:, k), :);
       end
       nrm = null(M);
       nrm = eps * nrm/norm(nrm, 2);
       
       %Get the point in both dirs then figure out which is extrenal
       new_p = avg_v' + nrm*[1, -1];
       %new_p = new_p(:, ~all(PolySet(i).H(:, 1:end-1) * new_p < PolySet(i).H(:, end)));
       pts(:, cur_pt_idx) = new_p(:, ~all(PolySet(i).H(:, 1:end-1) * new_p < PolySet(i).H(:, end)));
       cur_pt_idx = cur_pt_idx + 1;
    end
end

for i = 1:num_edges
    for j = 1:numel(PolySet)
        if(all(PolySet(j).H(:, 1:end-1) * pts(:, i) < PolySet(j).H(:, end)))
            pts(:, i) = nan * pts(:, i);
        end
    end
end
pts = pts(:, ~isnan(pts(1,:)));
end

function [k, s] = getLaw(pts, inputs)
K = inputs/[pts; ones(1, size(pts,2))];
k = K(:, 1:end-1);
s = K(:, end);
end