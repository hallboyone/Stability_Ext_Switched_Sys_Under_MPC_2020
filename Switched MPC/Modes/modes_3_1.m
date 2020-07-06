%================================ DYNAMICS ================================
numModes = 2;
m = [];
for i=1:numModes
    m = [m,sysMode];
end
%m = repmat(sysMode, 2, 1);

m(1).N = 5;
m(2).N = 4;

Ts1 = 0.1;
Ts2 = 0.15;
% assume two modes share the same state-space equation 
m(1).A = [1, Ts1, -Ts1; 0, 1, 0;0, 0, 1];
m(2).A = [1, Ts2, -Ts2; 0, 1, 0;0, 0, 1];

m(1).B = [Ts1^2/2, -Ts1^2/2; Ts1, 0; 0, Ts1];
m(2).B = [Ts2^2/2, -Ts2^2/2; Ts2, 0; 0, Ts2];

%============================== CONSTRAINTS ===============================
%                        Mode 1,      Mode 2,...
inputLowerBounds =   {[-0.3;-0.3], [-0.3;-0.3]};
inputUpperBounds =   {[0.1; 0.1], [0.1; 0.1]};

%The value of [MODE i L.B.] should be set to the 3x1 column vector where
%each element is the minimum value in that diminsion. 
%EXAMPLE
%stateLowerBounds = {-[20;10;2],-[10;5;1]};
%stateUpperBounds = { [20;10;2], [10;5;1]};

stateLowerBounds =  {[-0.5333; -0.2; -0.2], [-0.4; -0.2; -0.2]};
stateUpperBounds =  {[10; 0.6; 0.6], [5; 0.5; 0.5]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
%assume two modes share the same Q, R
m(1).Q = [0.5,0,0;0, 0.25,0;0,0,0.15];
m(1).R = 1.1*[1 0; 0 .75];

m(2).Q = [0.5,0,0;0, 0.25,0;0,0,0.5];
m(2).R = 1.1*[.75 0; 0 1];

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);