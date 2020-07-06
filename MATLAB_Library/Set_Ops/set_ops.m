classdef set_ops
    methods(Static)
        function pre_set = pre(mat, poly, steps)
            if ~exist('steps', 'var')
                steps = 1;
            end
            if size(mat,2)==1 && size(poly,2)==1 %Auto case
                A = mat{1};
                pre_set = poly{1};
                if ~set_ops.issquarematrix(A)
                    error("A must be a real square matrix");
                elseif ~isa(pre_set, 'Polyhedron')
                    error("X must be a Polyhedron");
                elseif pre_set.Dim~=size(A)
                    error("The diminsion of A and X must match");
                end
                pre_set = pre_set*A^steps;
            elseif size(mat,2)==2 && size(poly, 2)>=2 %Controlled case
                A = mat{1};
                B = mat{2};
                pre_set = poly{1};
                U = poly{2};
                
                if ~set_ops.issquarematrix(A)
                    error("A must be a real square matrix");
                elseif ~ismatrix(B) || ~isreal(B)
                    error("B must be a real matrix");
                elseif ~isa(pre_set, 'Polyhedron')
                    error("X must be a Polyhedron");
                elseif ~isa(U, 'Polyhedron')
                    error("U must be a Polyhedron");
                elseif pre_set.Dim~=size(A)
                    error("The diminsion of A and X must match");
                elseif U.Dim~=size(B, 2)
                    error("The diminsion U must match the number of columns in B");
                elseif size(A,1)~=size(B,1)
                    error("A and B must have the same number of rows");
                end
                
                BU = (-B)*U;
                BU.minVRep;
                BU.minHRep;
                for i=1:steps
                    temp_set = pre_set+BU;
                    temp_set.minHRep;
                    %temp_set.computeVRep;
                    pre_set = temp_set*A;
                    if(size(poly, 2) == 3)
                        pre_set = pre_set & poly{3};
                    end
                    %pre_set.computeVRep;
                end                    
            end
        end
        
        function reach_set = reach(mat, poly, steps)
           %mat - {A, opt-B}, poly - {targ, opt-U}}
           if ~exist('steps', 'var')
               steps = 1;
           end
           if size(mat,2)==1 && size(poly,2)==1 %Auto case
               A = mat{1};
               reach_set = poly{1};
               if ~set_ops.issquarematrix(A)
                   error("A must be a real square matrix");
               elseif ~isa(reach_set, 'Polyhedron')
                   error("The target set must be a Polyhedron");
               elseif reach_set.Dim~=size(A)
                   error("The diminsion of A and the target must match");
               end
               for i=1:steps
                   reach_set = A*reach_set;
               end
               
            elseif size(mat,2)==2 && size(poly, 2)==2 %Controlled case
                A = mat{1};
                B = mat{2};
                reach_set = poly{1};
                U = poly{2};
                if ~set_ops.issquarematrix(A)
                    error("A must be a real square matrix");
                elseif ~ismatrix(B) || ~isreal(B)
                    error("B must be a real matrix");
                elseif ~isa(reach_set, 'Polyhedron')
                    error("The target must be a Polyhedron");
                elseif ~isa(U, 'Polyhedron')
                    error("U must be a Polyhedron");
                elseif reach_set.Dim~=size(A)
                    error("The diminsion of A and the target must match");
                elseif U.Dim~=size(B, 2)
                    error("The diminsion U must match the number of columns in B");
                elseif size(A,1)~=size(B,1)
                    error("A and B must have the same number of rows");
                end
                for i=1:steps
                    reach_set = A*reach_set+B*U;
                    if rem(i,5)==0
                        reach_set.minHRep;
                    end
                end
           end
           
        end
        
        function cur_set = invar(mat, poly, starting_set, itr, verbose)
            if ~exist('itr', 'var')
                itr = 60;
            end
            if ~exist('starting_set', 'var')
                cur_set = poly{1}; %Set the starting set to all valid states
            else
                cur_set = starting_set; %Use the given starting set
            end
            if ~exist('verbose', 'var')
                verbose = 2;
            end
            
            if verbose > 1
                figure;
                hold on
            end
            for i=1:itr
                next_set = set_ops.pre(mat,{cur_set, poly{2:end}});
                next_set = next_set & cur_set;
                if verbose > 1
                    plot(cur_set, 'alpha', 0.1);
                    drawnow
                end
                if rem(i,3)==0 %Every fifth iteration, compute the minHRep
                    next_set = next_set.minHRep;
                end
                
                if next_set == cur_set
                    if verbose > 0
                        fprintf("Converged in %d iterations\n", i);
                    end
                    return;
                end
                cur_set = next_set;
            end
            warning("Did not converge");          
        end
        
        function result = issquarematrix(A)
            result = false;
            if ~ismatrix(A)     %If A is not a matrix, return false
                return;
            end
            
            if diff(size(A))~=0 %If A is not square, return false
                return;
            end
            
            if ~isreal(A)       %If A is not real, return false
                return;
            end
            
            result = true;
        end
    end
    
end