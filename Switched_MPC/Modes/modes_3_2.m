
%================================ DYNAMICS ================================
numModes = 2;
m = [];
for i=1:numModes
    m = [m,sysMode];
end

m(1).N = 10;
m(2).N = 10;

m(2).A = [1.913 0.06444 -0.002508; 
           41.24  1.913  -0.08357; 
           0, 0, 0.01111];
       
m(1).A = [2.496    0.07307  -0.001495; 
          71.61      2.496   -0.05506; 
          0          0   0.006738];

m(2).B = [-0.003613;-0.1881;0.8241];
m(1).B = [-0.00278;-0.1495;0.9933];

%============================== CONSTRAINTS ===============================
%                    Mode 1, Mode 2,...
inputLowerBounds =   {-0.75,-0.5};
inputUpperBounds =   { 0.75, 0.5};
stateLowerBounds = {-[20;20;20],-[20;20;20], -[20;20;20]};
stateUpperBounds = { [20;20;20], [20;20;20],  [20;20;20]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
m(1).Q = 0.1 * eye(3);
m(1).R = 0.01;

m(2).Q = 0.05 * eye(3);
m(2).R = 0.05;

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);



