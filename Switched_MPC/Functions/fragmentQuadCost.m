%Takes a sysMode object and polytope P and computes a cell array of cost matricies 
%which are valid within the repsective fragment of P
function [costs, fragments, fragK] = fragmentQuadCost(m, P)
[fragments, vSet] = makeWedges(P);
V = P.V;

dim = m(1).nx;
%Set of cost wieghts
costs = {repmat({nan(dim)}, numel(fragments), 1)};

%Inputs at each vertex
inputs = zeros(m.nu, m.N, size(V,1));
for i=1:size(V, 1)
    inputs(:,:,i) = cell2mat(m.MPC_Cntrl(V(i,:)'));
end

batchA = BatchA(m.A, m.N);
batchB = BatchB(m.A, m.B, m.N);
[batchQ, batchR] = BatchCost(m.Q, m.R, m.P, m.N);

%Cost mat M = K'AK + BK + C
A = batchB'*batchQ*batchB + batchR;
B = 2*batchA'*batchQ*batchB;
C = batchA'*batchQ*batchA;

Kset = zeros(m.nu*m.N, dim);
fragK = repmat({zeros(m.nu, dim)}, numel(fragments), 1);
nu = m.nu;
N = m.N;

for i=1:numel(fragments)
    F = inv(V(vSet(i, :), :)');
    for n=1:N
       U = reshape(inputs(:,n,vSet(i, :)),nu, dim);
       Kset((n-1)*nu+1:n*nu,:) = U*F;
       if n==1
           fragK{i} = U*F;
       end
    end
    costs{i} = Kset'*A*Kset + B*Kset + C;
end
end