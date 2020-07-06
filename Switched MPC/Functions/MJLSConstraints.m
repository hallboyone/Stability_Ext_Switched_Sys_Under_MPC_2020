function MJLSConstraints(sCost, gamma, name)
% MJLSCONSTRAINTS  Computes the probability and ADT given the switch costs and 
% dwell step decrease.

disp("The MJLS Constraints using the " + name + " method are:");

%Start with the lower and upper bounds of the probability
a = 1;
b = 0;

%Create the vector of dwell step decreases and switch costs
C = [gamma, sCost(end:-1:1)];
C((1+numel(gamma)):numel(gamma)+1:end) = [];

%Run a binary algorithm to up to 2^-16 accuracy
for i = 1:16
    %Get the value between the two bounds and build the standard P matrix
    p = (a+b)/2;
    P = p * eye(numel(gamma));
    P(P==0) = (1-p)/(numel(gamma)-1);
    
    %Create the transformed P matrix as laid out in paper
    P_transformed = transform(P);
    
    %Get the index of the eigenvalue equal to 1
    e_i = (abs(eig(P_transformed)-1) < 10e-10);
    
    %Get the cooresponding eigen vector
    [v, ~] = eig(P_transformed);
    v = v(:, e_i);
    
    %Normalize
    v = v/sum(v);
    
    %Check what side of the bound we are on and adjust accourdingly
    if C*v<1
        a = p;
    else
        b = p;
    end
end

%Print results
fprintf("Dwell probability - %0.6f\nAverage Dwell Time - %0.1f\n", p, p/(1-p));
end

function P = transform(P_orig)
%Transformed P matrix is comprized of 4 parts
% [[Diag of orig mat ], [Post Switch Region];
%  [Odds of switching], [   zeros region   ]]

%Get the number of modes
sz = size(P_orig);

%Create the dwell region (diag of orig mat)
dwell_region = P_orig;
dwell_region(eye(sz(1))==0) = 0;

%Create the post switch region
post_switch_region = repmat(eye(sz(1)), 1, sz(1));
for i = sz(1):-1:1
    post_switch_region(:,(i-1) * sz(1) + i) = [];
end

%Create the switch region building block
switch_region_block = P_orig;
switch_region_block(eye(sz(1)) == 1) = [];
switch_region_block = reshape(switch_region_block, sz(1) - 1, sz(1));

%Build the switch region
switch_region = repmat(switch_region_block, sz(1), 1);
for i = 1:sz(1)
    switch_region((i-1)*(sz(1)-1)+1:(i)*(sz(1)-1), [1:i-1,i+1:end]) = 0;
end

%Construct the full matrix
P = [dwell_region, post_switch_region; switch_region, zeros(size(switch_region, 1), size(post_switch_region, 2))];
end
