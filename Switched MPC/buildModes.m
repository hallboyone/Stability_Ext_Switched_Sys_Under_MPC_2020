function [m, costs] = buildModes(ops)

%% OPEN SYSTEM DEFINITION
%System diminsion and mesh shorthand
dim = ops.dim;
mesh = ops.mesh;

%Get an array of the indices of each point in the dim-diminsion grid
numPoints = mesh.numPoints;

%Make matrix of all X_vals at once. Greatly increases speed
X_vals = mesh.value(1:numPoints);

%Dummy var to silence warnings
m = 0;

if(ops.run_timer)
    tic;
end
%Load the basic mode data as sysMode objects
modeSet = str2func("modes_"+dim+"_"+ops.version);
modeSet();

%% LOAD SAVED MODE DATA
%Flags for if the mode needs to be rebuilt (i.e. if something changed)
needsFullRebuild = false(1,numModes);

%If a saved version exists, compare them to the new modes
if ~ops.run_timer && isfile("./Data/modes_"+dim+"_"+ops.version+".mat") && ~ops.force_rebuild_modes
    mSaved = load("./Data/modes_"+dim+"_"+ops.version+".mat");
    mSaved = mSaved.m(:);
    
    %If the number of modes don't match or we force rebuild
    if(numel(mSaved)~=numModes)
        needsFullRebuild(:) = true;
    else
        %Check each mode. If they are equivalent, get the saved set
        %results. Else, flag mode as needing set operations
        for i=1:numModes
            if m(i).equivalent(mSaved(i))
                m(i) = mSaved(i);
            else
                needsFullRebuild(i) = true;
            end
        end
    end
else
    needsFullRebuild(:) = true;
end

%% SET OPS & MPC SETUP
if(any(needsFullRebuild))
    %Get the terminal region and the set of pre sets for each mode which needs
    %to be rebuilt
    disp("Getting terminal sets...");
    for i=1:numModes
        if needsFullRebuild(i)
            m(i).setT();
        end
    end
    
    disp("Getting pre sets...");
    for i=1:numModes
        if needsFullRebuild(i)
            m(i).setPre();
        end
    end
    
    disp("Building MPC Controllers...");
    for i=1:numModes
        if needsFullRebuild(i)
            m(i).buildMPC(ops.Solver_ops);
        end
    end
end

%% BUILD COST STRUCTURE
%Build skeleton struct
costs = struct('Low', {repmat({[]}, numModes, 1)}, 'Norm', [], 'Valid', {repmat({[]}, numModes, 1)});

if(ops.run_timer)
    costs.time = 0;
end

%Add the field names the user wants to work with to the struct
fieldNames = fieldnames(ops);
for i=1:numel(fieldNames)
    if size(fieldNames{i},2)>8 && all(fieldNames{i}(1:9)=='Run_Cost_') && ops.(fieldNames{i})
        costs.(fieldNames{i}(10:end)) = repmat({[]}, numModes, 1);
    end
end
getTrueCost = isfield(costs, 'True');

%% LOAD AVALIBLE COST DATA
loadedValidCosts = false;
if ~ops.run_timer && isfile("./Data/costs_"+dim+"_"+ops.version+".mat")
    savedCosts = load("./Data/costs_"+dim+"_"+ops.version+".mat");
    savedMesh = savedCosts.mesh;
    savedCosts = savedCosts.costs;
    
    %Check if the mesh that was saved is the mesh that we are using
    if mesh.equals(savedMesh)
        %Check if the number of modes has changed
        if numel(savedCosts.Low) == numModes
            %If both are true, saved costs are valid
            loadedValidCosts = true;
            
            %Load the saved norm
            costs.Norm = savedCosts.Norm;
            
            %Load the save lower bound and feasible points
            for i=1:numModes
                costs.Low{i} = savedCosts.Low{i};
                costs.Valid{i} = savedCosts.Valid{i};
            end
            
            %Load the true cost if it was saved
            if isfield(savedCosts, 'True')
                for i=1:numModes
                    costs.True{i} = savedCosts.True{i};
                    costs.Inputs{i} = savedCosts.Inputs{i};
                end
            end
        end
    end
end
%Flags if a mode needs its costs recomputed. All true if no valid costs
%loaded. Else only changed modes are true
needsCostRebuild = needsFullRebuild | (~loadedValidCosts);

%% COMPUTE MISSING COST DATA
%Norms
if ~loadedValidCosts
    disp("Computing the norms...");
    costs.Norm = vecnorm(X_vals)';
end

%Lower Bound
if any(needsCostRebuild)
    %Get the feasible points. Points are stored as an array of their
    %indecies ([123, 124, 156, 157, 158,...])
    disp("Computing feasible points")
    feas_region = m(1).S(end);
    for i = 2:numModes
        feas_region = feas_region & m(1).S(end);
    end
    feas_region.minHRep;
    H = feas_region.H;
    for i=1:numModes
        %Get the H matrix defineing the feasible region
        %H = m(i).S(end).H;
        costs.Valid{i} = zeros(1, numPoints);
        
        %Look at the next 'step' points
        step = 500;
        cur_idx = 1:step;
        for point = 1:step:numPoints-step
            %For the group of points, get the ones which are feasible
            validPoints = cur_idx(all(H*[X_vals(:, point:point+step-1); -ones(1,step)] <=0));
            
            %Save the feasible points in their respective bins
            costs.Valid{i}(validPoints) = validPoints;
            cur_idx = cur_idx + step;
        end
        %Take care of the remaining points
        cur_idx = cur_idx(cur_idx<=numPoints);
        validPoints = cur_idx(all(H*[X_vals(:, cur_idx); -ones(1,numel(cur_idx))] <=0));
        costs.Valid{i}(validPoints) = validPoints;
        
        %Remove all 0 elements
        costs.Valid{i}(costs.Valid{i}==0) = [];
    end
    
    %Get the lower bound on the costs
    disp("Computing the lower bound...");
    for i=1:numModes
        if needsCostRebuild(i)
            %Initialise over entire mesh
            costs.Low{i} = nan(numPoints, 1);
            
            %Loop over just the valid points and save their LQR cost
            P = m(i).P;
            for idx = costs.Valid{i}
                costs.Low{i}(idx) = X_vals(:,idx)' * P * X_vals(:,idx);
            end
        end
    end
end

%True Cost
if getTrueCost
    for mode = 1:numModes
        if needsCostRebuild(mode)
            fprintf("Computing the true cost of mode %d...\n", mode);
            costs.True{mode} = nan(numPoints, 1);
            costs.Inputs{mode} = nan( m(mode).nu, numPoints);
            for i = costs.Valid{mode}
                %x = X_vals(:, i);
                %Compute the optimal cost
                [costs.True{mode}(i), costs.Inputs{mode}(:, i)] = m(mode).MPCCost(X_vals(:, i));%quadCost(@(x,u) m(mode).f(x,u), x, cell2mat(m(mode).MPC_Cntrl(x)),...
                    %m(mode).Q, m(mode).R, m(mode).P);
            end
        end
        %Set the cost at the origin to 0
        try
            costs.True{mode}(mesh.origin) = 0;
            costs.Inputs{mode}(mesh.origin) = costs.Inputs{mode}(mesh.origin) * 0;
        catch
            warning("Without an origin, there may be unexpected behavior")
        end
    end
end

costs.Norm = round(costs.Norm, 10);
if(ops.run_timer)
    costs.time = toc;
end

%% SAVE NEW DATA
save("./Data/modes_"+dim+"_"+ops.version, "m");
save("./Data/costs_"+dim+"_"+ops.version, 'costs', 'mesh');

fprintf("Loaded v.%d of %dD modes\n", ops.version, dim);
end







