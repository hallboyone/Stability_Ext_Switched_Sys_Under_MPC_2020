
%================================ DYNAMICS ================================
numModes = 3;
m = [];
for i=1:numModes
    m = [m,sysMode];
end

m(1).N = 6;
m(2).N = 8;
m(3).N = 10;

m(1).A = [1, .1; -.1, 1.];
m(2).A = [1, -.1; .05, 1];
m(3).A = [1, 0.2; 0 1];

m(1).B = [0.0;.25];
m(2).B = [0.0;.5];
m(3).B = [0; .15];

%============================== CONSTRAINTS ===============================
%                    Mode 1, Mode 2,...
inputLowerBounds =   {-0.75,-0.5, -1};
inputUpperBounds =   { 0.75, 0.5, 1};
stateLowerBounds = {-[5;5],-[5;5], -[5;5]};
stateUpperBounds = { [5;5], [5;5], [5;5]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
m(1).Q = [1 0; 0 5];
m(1).R = 1;

m(2).Q = [5 0; 0 1];
m(2).R = 1;

m(3).Q = eye(2);
m(3).R = 1;

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);