function Unstable_3D_2()

A = {};
B = {};
X = {};
U = {};
Q = {};
R = {};

Ts1 = 0.1;
Ts2 = 0.15;

% equ_points = {[3; 1; 1], [2; 0.5; 0.5]};

% assume two modes share the same state-space equation 
A{1} = [1, Ts1, -Ts1; 0, 1, 0;0, 0, 1];
A{2} = [1, Ts2, -Ts2; 0, 1, 0;0, 0, 1];

B{1} = [Ts1^2/2, -Ts1^2/2; Ts1, 0; 0, Ts1];
B{2} = [Ts2^2/2, -Ts2^2/2; Ts2, 0; 0, Ts2];


inputLowerBounds =   {[-0.3;-0.3], [-0.3;-0.3]};
inputUpperBounds =   {[0.1; 0.1], [0.25; 0.25]};


stateLowerBounds =  {[-0.5333; -0.4; -0.4], [-0.4; -0.2; -0.2]}; 
stateUpperBounds =  {[10; 2; 2], [5; 1; 1]}; 

Q = {};
R = {};
lim = [0, 1.23];

Q{1} = [2 -1 0; -1 2 -1;0 1 2];
R{1} = [1 0; 0 .75];

Q{2} = [0.3,0,0;0, 0.6,0;0,0,0.9];
R{2} = [0.75 0; 0 1];

for i=1:2
%     [Q{i}, R{i}] = cost_function(3, 2, lim(i));
    X{i} = Polyhedron([eye(3); -eye(3)],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    U{i} = Polyhedron([eye(2); -eye(2)],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

X_feasible = X{1}.intersect(X{2});

N = 5;

controllers = {};
for i = 1:2
    controllers{end+1} = buildcontrollers(A{i}, B{i}, X{i}, U{i}, Q{i}, R{i}, N);
end

%%
% [mu, s, J1, J2] = simulations(controllers, A, B, Q, R, X_feasible, equ_points, 100*N);
sim_range = 2000*N;

mu = [];
J1 = [];
J2 = [];

mode = 2; %begin with mode 2
numStep = N; 

% generate a random initial point which is included by X feasible set
x = [0.831817912727996;0.823204476836442;0.385800245714414]; %-> like this initial point will be unstable
% while numStep > 0 
%     x = (3).*rand(3,1);
%     if X_feasible.isInside(x) == 1
%         break;
%     end
% end

T = 0;
s = 1;
mu = [mu; (x)'];
J1 = [J1, calculate_cost(Q{1}, R{1}, x, controllers{1})];



for n = 1:sim_range
    numStep = numStep - 1;
    u = controllers{mode}{x};
   
    if isnan(u)
        disp("no u, infeasible")
        break;
    end
      
%     if numStep < 0
%         if mode == 1
%             J1 = [J1, (x)'*Q{mode}*(x) + u'* R{mode}*u];
%             u2 = controllers{2}{x};
%             if ~isnan(u2)
%                 [mode, J] = check_switch(Q, R, mode, x, u2, J2);
%                 if mode ~= 2
%                     J2 = J;
%                     numStep = N;
%                 end
%             end
%         else
%             J2 = [J2, (x)'*Q{mode}*(x) + u'* R{mode}*u];
%             u1 = controllers{1}{x};
%             if ~isnan(u1)
%                 [mode, J] = check_switch(Q, R, mode, x, u1, J1);
%                 if mode ~= 1
%                     J1 = J;
%                     numStep = N;
%                 end
%             end
%         end
%     end
         
    x = A{mode} * x + B{mode} * u; 
          
    mu = [mu; (x)'];
        
end
if n == sim_range
    s = 0;
end


%%
figure(1)
subplot(3, 2, [1, 3, 5])
plot3(mu(:, 1), mu(:, 2), mu(:, 3))
grid on 
xlabel('d')
ylabel('v1')
zlabel('v2')
subplot(3, 2, 2)
plot(mu(:, 1))
ylabel('d')
subplot(3, 2, 4)
plot(mu(:, 2))
ylabel('v1')
subplot(3, 2, 6)
plot(mu(:, 3))
ylabel('v2')

end


function controllers = buildcontrollers(A, B, X, U, Q, R, N)

ops = sdpsettings('verbose', 0, 'CACHESOLVERS', 1);


x = sdpvar(repmat(3, 1, N+1), ones(1, N+1)); %State variable
u = sdpvar(repmat(2, 1, N), ones(1, N));     %Input variable

objective = x{N+1}'*Q*x{N+1};
constraints = [];

for k=1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1}==A*x{k} + B*u{k},...     %Dynamic Constraint
        X.A*x{k}<=X.b,...%State Constraint
        U.A*u{k}<=U.b];  %Input Constraint
end

controllers = optimizer(constraints, objective ,ops,x{1},u{1});

end


function J = calculate_cost(Q, R, x, controller)
    u = controller{x}; 
    J = x'*Q*x + u'*R*u;
end

function [mode, J] = check_switch(Q, R, curr_mode, x, u, J)

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
    if (x'*Q{2}*x + u'* R{2}*u) > max(J)
        mode = 2;
        disp("change to mode 2")
        J = [J, x'*Q{2}*x + u'* R{2}*u];
    else
        mode = 1;
    end
else
    if (x'*Q{1}*x + u'* R{1}*u) > max(J)
        mode = 1;
        disp("change to mode 1")
        J = [J, x'*Q{1}*x + u'* R{1}*u];
    else
        mode = 2;
    end
end
  
end




