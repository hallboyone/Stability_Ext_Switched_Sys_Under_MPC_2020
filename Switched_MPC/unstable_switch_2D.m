function unstable_switch_2D()

% 2D case 
A = {};
B = {};
X = {};
U = {};
Q = {};
R = {};

% assume two modes share the same state-space equation 
A{1} = [-0.1, -1; 2, -0.1];
A{2} = [-0.1, 2;-1, -0.1];

B{1} = [0.5, 0.1];
B{2} = [-1, 3];


inputLowerBounds =   {[-0.3;-0.3], [-0.3;-0.3]};
inputUpperBounds =   {[0.1; 0.1], [0.1; 0.1]};


stateLowerBounds =  {[-20; -20], [-30; -10]}; 
stateUpperBounds =  {[20; 20], [10; 30]}; 

Q{1} = [0.3,0;0, 0.9];
R{1} = [1 0; 0, 0.5];

Q{2} = [0.9, 0; 0, 0.3];
R{2} = [0.5 0; 0 1];

for i=1:2
    X{i} = Polyhedron([eye(2); -eye(2)],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    U{i} = Polyhedron([eye(2); -eye(2)],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

X_feasible = X{1}.intersect(X{2});

N = 10;
controllers = buildcontrollers(A, B, X, U, Q, R, N);

mu = simulations(controllers, A, B, Q, R, X_feasible, 1000*N);

plot3(mu(:, 1), mu(:, 2), zeros(length(mu(:, 1)), 1))

end


function controllers = buildcontrollers(A, B, X, U, Q, R, N)

controllers = {};

ops = sdpsettings('verbose', 0, 'CACHESOLVERS', 1);

for i = 1:2
    x = sdpvar(repmat(2, 1, N+1), ones(1, N+1)); %State variable
    u = sdpvar(repmat(2, 1, N), ones(1, N));     %Input variable
    
    objective = x{N+1}'*Q{i}*x{N+1};
    constraints = [];
    
    for k=1:N
        objective = objective + x{k}'*Q{i}*x{k}+ u{k}'*R{i}*u{k};
        constraints = [constraints, x{k+1}==A{i}*x{k} + B{i}*u{k},...     %Dynamic Constraint
            X{i}.A*x{k}<=X{i}.b,...%State Constraint
            U{i}.A*u{k}<=U{i}.b];  %Input Constraint
    end      
    controllers{i} = optimizer(constraints, objective ,ops,x{1},u{1});
    
end

end


function mu = simulations(controllers, A, B, Q, R, X, sim_range)

mu = [];
J1 = [];
J2 = [];

step = sim_range/10;

mode = 2; %begin with mode 2

% generate a random initial point which is included by X feasible set
while step > 0 
    x = (0.01).*rand(2,1);
    if X.isInside(x) == 1
        break;
    end
end

mu = [mu; x'];

for n = 1:sim_range
    u = controllers{mode}{x};
   
    if isnan(u)
        disp("1, no u, infeasible")
        break;
    end
    
    
    if mode == 1
        u2 = controllers{2}{x};
        if isnan(u2)
            disp("2, no u, infeasible")
            clc
            break;
        else
            mode = check_switch(Q, R, mode, x, u2, J2);
            if mode ~= 1
               disp("switch to mode 2")
            end
            J2 = [J2, x'*Q{2}*x + u2'* R{2}*u2];
        end
        J1 = [J1, x'*Q{1}*x + u'* R{1}*u];
    else
        u1 = controllers{1}{x};
        if isnan(u1)
            disp("3, no u, infeasible")
            clc
            break;
        else
            mode = check_switch(Q, R, mode, x, u1, J1);
            if mode ~= 2
                disp("switch to mode 1")
            end
            J1 = [J1, x'*Q{1}*x + u1'* R{1}*u1];
        end 
        J2 = [J2, x'*Q{2}*x + u'* R{2}*u];
    end
    
    x = A{mode} * x + B{mode} * u; 
    
       
    mu = [mu; x'];      
end

end

function mode = check_switch(Q, R, curr_mode, x, u, J)

% attempt to use Lypunov equation to check 
% Basic proccess: satrt with mode 2 and x1
% calculate u1 and u2 in mode 1 controller and mode 2 controller
% record the first x1 energy (x'Qx + u'Ru) in J1 and J2.
% use A x1+ B u1 = x to next state
% calculate x energy in mode 1 and 2
% if x energy in mode 1 is larger than the all J1 energy, switch to mode 1
% else continue. 
% repeat 

if curr_mode == 1
    J_2 = x'*Q{2}*x + u'* R{2}*u;
    if J_2 > max(J)
        mode = 2;
        disp("change to mode 2")
    else
        mode = 1;
    end
else
    J_1 = x'*Q{1}*x + u'* R{1}*u;
    if J_1 > max(J)
        mode = 1;
        disp("change to mode 1")
    else
        mode = 2;
    end
end
  
end








