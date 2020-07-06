
%================================ DYNAMICS ================================
numModes = 2;
m = [];
for i=1:numModes
    m = [m,sysMode];
end

m(1).N = 5;
m(2).N = 5;

m(1).A = [0.8,-1.0;
          1.0, 0.7];
      
m(2).A = [1.2, 0.5;
          0.0, 1.2];

m(1).B = [0; 1.0];
m(2).B = [0; 0.5];

%============================== CONSTRAINTS ===============================
%                    Mode 1, Mode 2,...
inputLowerBounds =   {-340, -600};
inputUpperBounds =   { 340,  600};
stateLowerBounds = {-2000*[1;1],-1600*[1;1]};
stateUpperBounds = { 2000*[1;1], 1600*[1;1]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
m(1).Q = 10 * eye(2);
m(1).R = 1;

m(2).Q = 5 * eye(2);
m(2).R = 1;

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);