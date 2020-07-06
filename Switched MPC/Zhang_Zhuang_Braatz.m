function Zhang_Zhuang_Braatz(c, modes, ops)
num_m = numel(modes);
run_cor = [0,0,1,1];

if (ops.run_timer)
    tic;
end
%Generate the intersection of all the feasible and terminal regions along
%with the dynamic matricies under LQR control
joint_feas_set = modes(1).X;
joint_term_set = modes(1).X;
A_hats         = cell(num_m, 1);
for m_i = 1:num_m
    joint_feas_set = modes(m_i).S(end) & joint_feas_set;
    joint_term_set = modes(m_i).T      & joint_term_set;
    A_hats{m_i} = modes(m_i).A - modes(m_i).B*modes(m_i).K;
end

tic
if (ops.dim==2 && ops.version==3)
    Delta = [2,1];
    O = Algorithm1A(joint_term_set, A_hats, Delta, 0.99, 0);
elseif (ops.dim==3 && ops.version==1)
    Delta = [1,1];
    O = Algorithm1A(joint_term_set, A_hats, Delta, 0.99, 0);
elseif (ops.dim==2 && ops.version==4)
    Delta = [2,1,2];
    O = Algorithm1A(joint_term_set, A_hats, Delta, 0.8, 0);
elseif (ops.dim==2 && ops.version==5)
    Delta = [1,1];
    O = Algorithm1A(joint_term_set, A_hats, [1,1], 0.9, 0);
end
O_time = toc;

%% Find the MDTs based on Cor 2
if (run_cor(2) || run_cor(1) || run_cor(3))
    fprintf("Running Cor 2...\n");
    MDT_Cor_2 = nan(numel(modes), 1);
    for m_i=1:num_m
        %Get the feasible region of the mode
        S = modes(m_i).S(end);
        
        %Init dwell time to 0
        MDT_Cor_2(m_i, 1) = 0;
        
        %While the mode has not dwelled long enough to ensure we are in the safe
        %region
        while(~joint_feas_set.contains(S))
            %Compute the 1-step reach set under MPC from current set
            S = modes(m_i).MPCReachSet(S, 1);
            %Increment dwell time
            MDT_Cor_2(m_i, 1) = MDT_Cor_2(m_i, 1) + 1;
        end
        S.minVRep();
    end
end
%% Compute the MDTs based on Cor 1
if (run_cor(1))
    MDT_Cor_1 = nan(numel(modes), 2);
    fprintf("Running Cor 1...\n");
    %The first-stage dwell times equal those from Cor 2
    MDT_Cor_1(:,1) = MDT_Cor_2;
    for m_i=1:num_m
        %Get an array of the reach sets from the target set
        reach_sets = repmat(Polyhedron, 1, MDT_Cor_1(m_i, 1)+1);
        reach_sets(1) = joint_feas_set;
        for i=1:numel(reach_sets)-1
            reach_sets(i+1) = modes(m_i).MPCReachSet(reach_sets(i), 1);
        end
        reach_sets(1) = [];
        
        %Find the LAST reach set not in the target set
        last_set_i = MDT_Cor_1(m_i, 1);
        while(joint_feas_set.contains(reach_sets(last_set_i)))
            last_set_i = last_set_i-1;
            if last_set_i == 0
                break;
            end
        end
        
        %The system must dwell till it is past the last reach set
        MDT_Cor_1(m_i, 2) = last_set_i + 1;
    end
end

if (run_cor(3))
    fprintf("Running Thm 3...\n");
    MDT_Cor_3 = ThmThree(modes, MDT_Cor_2, O);
    MDT_Cor_3(:,end) = Delta;
end
%% Compute the MDTs based on Cor 4
if (run_cor(4))
    fprintf("Running Cor 4...\n");
    MDT_Cor_4 = nan(num_m, 2);
    MDT_Cor_4(:, 2) = Delta;
    
    for m_i=1:num_m
        dwell_time = 0;
        S = modes(m_i).S(end);
        while(~O.contains(S))
            dwell_time = dwell_time + 1;
            S = modes(m_i).MPCReachSet(S, 1);
        end
        MDT_Cor_4(m_i, 1) = dwell_time;
    end
end

if(ops.run_timer)
    c.time = c.time + toc;
    c.time
end
%% Print results
if(run_cor(1))
    fprintf("Corollary 1\n");
    for i=1:num_m
        fprintf(" - Mode %d MDT: %d (Stage 1), %d (Stage 2+)\n", i, MDT_Cor_1(i, 1), MDT_Cor_1(i, 2));
    end
end

if(run_cor(2))
fprintf("\nCorollary 2\n");
for i=1:num_m
    fprintf(" - Mode %d MDT: %d\n", i, MDT_Cor_2(i));
end
end

if(run_cor(3))
    fprintf("\nCorollary 3\n");
    fprintf("Stage#|");
    for i=1:size(MDT_Cor_3, 2)
        if (i<10)
            fprintf("___%d", i);
        else
            fprintf("__%d", i);
        end
    end
    fprintf("_\n");
    for i=1:num_m
        fprintf("Mode %d|", i);
        for n=1:size(MDT_Cor_3, 2)
            fprintf("%4d", MDT_Cor_3(i,n));
        end
        fprintf("\n");
    end
end

if(run_cor(4))
    fprintf("\nCorollary 4\n");
    for i=1:num_m
        fprintf(" - Mode %d MDT: %d (Stage 1), %d (Stage 2+)\n", i, MDT_Cor_4(i, 1), MDT_Cor_4(i, 2));
    end
end
end

function MDT_Thm_3 = ThmThree(mode, tau_first_stage, O)
% tau_first_stage   - MDT st all feas i.c. enter the inter of feas regions
% tau_i, i in [2,v] - MDT st (14) holds
% tau_i, i > v      - MDT which generated DT-contractive set O

v = 2; % At what stage do we enter the contractive set O?
while(true)
    joint_set = mode(1).MPCReachSet(mode(1).S(end), v-1);
    for m_i=2:numel(mode)
        joint_set = joint_set & mode(m_i).MPCReachSet(mode(m_i).S(end), v-1);
    end
    if O.contains(joint_set)
        break;
    else
        v = v+1;
    end
end
fprintf("Thm 3 requires %d M-MDT sets\n", v);

MDT_Thm_3 = nan(numel(mode), v+1);
MDT_Thm_3(:,1) = tau_first_stage;

starting_set = mode(1).S(end);
for m_i = 2:numel(mode)
    starting_set = starting_set & mode(m_i).S(end);
end
starting_set.minHRep();

for l = 2:v
    l
    target_set   = mode(1).MPCReachSet(mode(1).S(end), l-1);
    for m_i = 2:numel(mode)
        target_set = target_set & mode(m_i).MPCReachSet(mode(m_i).S(end), l-1);
    end
    target_set.minHRep();
    
    for m_i = 1:numel(mode)
        %Get an array of the reach sets from the starting set
        reach_sets = repmat(Polyhedron, 1, tau_first_stage(m_i)+1);
        reach_sets(1) = starting_set;
        for i = 2:numel(reach_sets)
            reach_sets(i) = mode(m_i).MPCReachSet(reach_sets(i-1), 1);
        end
        
        %Find the LAST reach set not in the target set
        last_set_i = tau_first_stage(m_i)+1;
        while(target_set.contains(reach_sets(last_set_i)) || O.contains(reach_sets(last_set_i)))
            last_set_i = last_set_i-1;
        end
        
        %The system must dwell till it is past the last reach set
        MDT_Thm_3(m_i, l) = last_set_i;
    end
    
    starting_set = target_set;
end
end

function O = Algorithm1A(X, A, Delta, gamma, show_plots)
% Implemanetation of Algorithm 1A in Zhang2015

% Get O_1
for m = 1:numel(A)
    S = preSetIntersection(A{m}, X, 0, Delta(m)-1);
    if (m==1)
        %NOTE: This is diffenent that in Zhang2015. They do not intersect
        %it with the terminal regions but this makes no sense. What that
        %outputs is an O-set that is not a subset of X. 
        O = X & preSetIntersection(A{m}, gamma*S, Delta(m), Delta(m)*2-1);
    else
        O = O & preSetIntersection(A{m}, gamma*S, Delta(m), Delta(m)*2-1);
    end
end
O.minHRep();

% Get the first iteration's result O_2
O_next = O;
for m = 1:numel(A)
    O_next = O_next & preSetIntersection(A{m}, gamma*O, Delta(m), 2*Delta(m)-1);
end
O_next.minHRep();

if show_plots
    figure;
    hold on;
    plot(O);
    drawnow;
end

% Untill the iterations converge, run the algorithm
while(O ~= O_next)
    O = O_next;
    for m = 1:numel(A)
        O_next = O_next & preSetIntersection(A{m}, gamma*O, Delta(m), 2*Delta(m)-1);
    end
    O_next.minHRep();
    
    if show_plots
        plot(O);
        drawnow
    end
end
end

function Y = preSetIntersection(A, X, t_low, t_high)
Y = X*A^t_low;
for i=(t_low+1):(t_high)
    S = X*(A^i);
    Y = Polyhedron([Y.H(:,1:end-1); S.H(:,1:end-1)], [Y.H(:,end); S.H(:,end)]); 
end
Y.minHRep();
end