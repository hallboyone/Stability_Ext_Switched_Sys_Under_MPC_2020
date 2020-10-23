
%% Building option struct
% Set Run_* to 0 or 1 to either skip or run the indicated script. The dim and
% version fields are used to select the correct script with the name
% modes_[dim]_[ver].m which loads the modes into the workspace. 
tic;
%Input program parameters. Indicate which UB to calculate
options = struct(...
    'Run_Cost_True',           1,...
    'Run_Cost_Fragments',      1,... 
    'Run_Sandbox',             0,...
    'Run_Zhang_Zhuagn_Braatz', 0,...
    'force_rebuild_modes',     0,... %Force everything to rebuild
    'dim',                     2,... %Diminsion of the modes
    'version',                 4,... %The version of modes in set diminsion
    'figs',                    0,... %Should figure be plotted.
    'Solver_ops',    sdpsettings(... %Options used by yalmip
                    'verbose', 0,...
               'solver','mosek'),...
    'DEBUG',                   0,... %Extra error checks
    'run_timer',               0);

options = checkFigs(options);
options = checkTimer(options);
options = setGridEdge(options);
options = makeMesh(options);

%Make structures for the gamma values of each mode and the S_{i,j} values
%for each switch
gamma  = struct;
S_cost = struct;

%% Load the model data
% This includes
%
% # Dynamics, constraints, and costs
% # MPC controllers
% # Terminal and pre sets
% # The norms, true cost, and lower bound at each grid point
% # The points at which mode i is feaible

[m, c] = buildModes(options);

%% Run Active Scripts
% Call various functions to run aleye(2)gorithms related to switched MPC on the
% modes using the program options

if(options.Run_Cost_True)
    [gamma.True, S_cost.True] =...
        TrueParameters(c, m, options);
end

if(options.Run_Cost_Fragments)
    [c.Fragments, gamma.Fragments, S_cost.Fragments] =...
        FragmentsUB(c, m, options);
end

if options.Run_Zhang_Zhuagn_Braatz
    Zhang_Zhuang_Braatz(c, m, options);
end

if options.Run_Sandbox
    sandbox(c, m, options);
end
toc
%% Helper Functions for Setup
% Helper functions used to set depended option values and check that values
% are valid. 

function ops = checkFigs(ops)
% CHECKFIGS Makes sure the figs field is only active when the diminsion of
% the system is less than 3
ops.figs = ops.figs && (ops.dim<=2);
end

function ops = setGridEdge(ops)
% SETGRIDEDGE Uses the version and dim feilds in ops to set the correct grid_edge
% field.
if (ops.dim==3 && ops.version == 1)
    ops.grid_edge = {-0.6:0.1:10,... %3.1
                     -0.2:.01:0.6,...
                     -0.2:.01:0.6};
elseif (ops.dim==2 && (ops.version == 1 || ops.version == 2))
   ops.grid_edge = -5:0.025:5;
elseif (ops.dim==2 && ops.version == 3)
    ops.grid_edge = -1300:10:1300;
elseif (ops.dim==2 && ops.version == 4)
    ops.grid_edge = -10:.025:10;
elseif (ops.dim==2 && ops.version == 5)
    ops.grid_edge = -50:.5:50;
end
end
%     'grid_edge',      -0.5:0.05:3,...

function ops = makeMesh(ops)
% MAKEMESH Uses the grid_edge field to construct a DiscreteGrid object to
% sample points at.
if numel(ops.grid_edge)==1 || ~iscell(ops.grid_edge)
    ops.mesh = DiscreteGrid(ops.grid_edge, ops.dim);
else
    ops.mesh = DiscreteGrid(ops.grid_edge);
    if ops.mesh.dim ~= ops.dim
        error("Diminsion of mesh is larger than diminsion of system!");
    end
end
end

function ops = checkTimer(ops)
% CHECKTIMER If the run_timer field is set and, if so, makes sure that only
% one function is called.
if (ops.run_timer)
    fieldNames = fieldnames(ops);
    num2run = 0;
    for i=1:numel(fieldNames)
        if size(fieldNames{i},2)>8 && all(fieldNames{i}(1:9)=='Run_Cost_') && ops.(fieldNames{i})
            if num2run == 0
                num2run = num2run + 1;
            else
                ops.(fieldNames{i}) = 0;
            end
        end
    end
end
end