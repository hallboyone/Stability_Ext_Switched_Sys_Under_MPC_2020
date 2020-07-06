
%================================ DYNAMICS ================================
numModes = 2;
m = [];
for i=1:numModes
    m = [m,sysMode];
end

m(1).N = 3;
m(2).N = 3;

m(1).A = [1, 1.5;
          0, 0.75];
      
m(2).A = [0.75, 0;
          1.5, 1];

m(1).B = [30; 0];
m(2).B = [0; 30];

%============================== CONSTRAINTS ===============================
%                    Mode 1, Mode 2,...
inputLowerBounds =   {-0.1, -0.1};
inputUpperBounds =   { 0.1,  0.1};
stateLowerBounds = {-50*[1;1],-50*[1;1]};
stateUpperBounds = { 50*[1;1], 50*[1;1]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
m(1).Q = [1 0; 0 10];
m(1).R = 1;

m(2).Q = [10 0; 0 1];
m(2).R = 1;

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);