classdef sysMode < matlab.mixin.Copyable
    properties
        A;
        B;
        Q;
        R;
        P;
        K;
        X;
        X_LQR;
        U;
        T;
        C; 
        S; %Pre sets
        N;
        
        MPC_Cntrl; %Returns all states
        MPC_u;     %Returns initial input only
        cur_x;
    end
    methods
        function x = f(obj, x, u)
            for i=1:size(u,2)
                x = obj.A*x + obj.B*u(:,i);
            end
        end
        
        function x = fMPC(obj, x)
            x = obj.A*x + obj.B*obj.MPC_u(x);
        end
        function isFeasible = inputFeasible(obj, x, u)
            for i=1:size(u,2)
                if ~obj.X.contains(x) || ~obj.U.contains(u(:,i))
                    isFeasible = false;
                    return;
                end
                x = obj.A*x + obj.B*u(:,i);
            end
            isFeasible = obj.T.contains(x);
        end
        
        function input_size = nu(obj)
            input_size = size(obj.B, 2);
        end
        
        function state_size = nx(obj)
            state_size = size(obj.A, 2);
        end
        
        function obj = idare(obj)
            [obj.P, obj.K, ~] = idare(obj.A, obj.B, obj.Q, obj.R);
            
        end
        
        function obj = setT(obj)
            obj.X_LQR = Polyhedron(obj.U.H(:,1:end-1) * (-obj.K), obj.U.H(:,end)) & obj.X;
            obj.T = set_ops.invar({obj.A-obj.B*obj.K}, {obj.X_LQR}, obj.X_LQR, 60, 0);
        end
        
        function obj = setC(obj)
            obj.C = set_ops.invar({obj.A, obj.B}, {obj.X, obj.U}, obj.X, 60, 0) & obj.X;
        end
        
        function obj = setPre(obj)
            obj.S = repmat(Polyhedron, obj.N+1, 1);
            obj.S(1) = obj.T;
            for i=1:obj.N
                obj.S(i+1) = set_ops.pre({obj.A, obj.B}, {obj.S(i), obj.U}) & obj.X;
                
                %For some reason, without this line the output of this
                %function was not defined.
                obj.S(i+1).minVRep;
                obj.S(i+1).minHRep;
            end
        end
        
        function obj = buildMPC(obj, ops)
            x = sdpvar(repmat(obj.nx, 1, obj.N+1), ones(1, obj.N+1)); %State variable
            u = sdpvar(repmat(obj.nu, 1, obj.N), ones(1, obj.N));     %Input variable
            
            objective = x{obj.N+1}'*obj.P*x{obj.N+1};
            constraints = obj.T.H(:, 1:end-1)*x{obj.N+1}<=obj.T.H(:, end);
            
            for k=1:obj.N
                objective = objective + x{k}'*obj.Q*x{k} + u{k}'*obj.R*u{k};
                constraints = [constraints, x{k+1}==obj.A*x{k} + obj.B*u{k},...     %Dynamic Constraint
                    obj.X.H(:, 1:end-1)*x{k}<=obj.X.H(:, end),...%State Constraint
                    obj.U.H(:, 1:end-1)*u{k}<=obj.U.H(:, end)];  %Input Constraint
            end
            
            %Optimizer object. First state given. Find optimal input.
            %Return all inputs
            obj.MPC_Cntrl = optimizer(constraints, objective ,ops,x{1},u);
            obj.MPC_u = optimizer(constraints, objective ,ops,x{1},u{1});
        end
        
         
        function [cost, input] = MPCCost(obj, x)
            %Returns the optimal cost of the system from state x
            u = cell2mat(obj.MPC_Cntrl(x));
            input = u(:,1);
            cost = quadCost(@(x,u) obj.f(x,u), x, u, obj.Q, obj.R, obj.P);
        end
        
        function R = MPCReachSet(obj, starting_set, step_count)
            %Returns the set that S maps to under MPC control
            
            V = starting_set.V'; %Make array of verticies
            
            for step_i=1:step_count
                for v_i = 1:size(V, 2)
                    if obj.T.contains(V(:,v_i))
                        V(:,v_i) = (obj.A - obj.B * obj.K) * V(:, v_i);
                    else
                        V(:,v_i) = obj.A*V(:,v_i) + obj.B * obj.MPC_u(V(:,v_i));
                    end
                end
            end
            
            V = V(:, any(~isnan(V)));
            R = Polyhedron(V');
            R.minVRep();
        end
        
        function equal = equivalent(obj, rhs)
            equal = true;
            if ~isequal(obj.A, rhs.A)
                equal = false;
            elseif ~isequal(obj.B, rhs.B)
                equal = false;
            elseif ~isequal(obj.Q, rhs.Q)
                equal = false;
            elseif ~isequal(obj.R, rhs.R)
                equal = false;
            elseif obj.X ~= rhs.X
                equal = false;
            elseif obj.U ~= rhs.U
                equal = false;
            elseif obj.N ~= rhs.N
                equal = false;
            end
        end
    end
end
    