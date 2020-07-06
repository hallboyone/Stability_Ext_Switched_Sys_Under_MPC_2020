function cost = crossCost(dyn_fun, x_1, x_2, u_1, u_2, Q, R, P)
% CROSSCOST  Computes the quadratic cost of the lin 
%   C = mpcCost(f, x, u, Q, R, P) simulates the cost of the sequence f(x,u)
%   using the cost function:
%   J(x, u) = x_n'*P*x_n + sum_(i=0)^(n-1) [x_i'*Q*x_i + u_i'*R*u_i]

if(isrow(inputs))
    inputs = inputs';
end

cost = 0;
for i=1:size(inputs, 1)
    cost = cost + cur_x'*Q*cur_x + inputs(i)'*R*inputs(i);
    cur_x = dyn_fun(cur_x, inputs(i));
end
cost = cost + cur_x'*P*cur_x;
end