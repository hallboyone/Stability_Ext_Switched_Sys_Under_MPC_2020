function cost = quadCost(dyn_fun, cur_x, inputs, Q, R, P)
% QUADCOST  Computes the quadratic cost of a series of inputs given the
% starting location 
%   C = mpcCost(f, x, u, Q, R, P) simulates the cost of the sequence f(x,u)
%   using the cost function:
%   J(x, u) = x_n'*P*x_n + sum_(i=0)^(n-1) [x_i'*Q*x_i + u_i'*R*u_i]

if(size(inputs, 1) ~= size(R, 1))
    inputs = inputs';
end

cost = 0;
for i=1:size(inputs, 2)
    cost = cost + cur_x'*Q*cur_x + inputs(:,i)'*R*inputs(:,i);
    cur_x = dyn_fun(cur_x, inputs(:, i));
end
cost = cost + cur_x'*P*cur_x;
end