%================================ DYNAMICS ================================
numModes = 3;
m = [];
for i=1:numModes
    m = [m,sysMode];
end

m(1).N = 4;
m(2).N = 4;
m(3).N = 4;

m(1).A = [0.0, 1.0;
         -2.5, 3.2];
      
m(2).A = [0.0, 1.0;
         -4.3, 4.5];
      
m(3).A = [0.0,  1.0;
          5.3, -5.2];

m(1).B = [0; 1];
m(2).B = [0; 1];
m(3).B = [0; 1];

%============================== CONSTRAINTS ===============================
%                    Mode 1, Mode 2,...
inputLowerBounds =   {-2.5, -5, -3.5};
inputUpperBounds =   { 2.5,  5,  3.5};
stateLowerBounds = {-30*[1;1],5*[-8;-10], 5*[-12; -10]};
stateUpperBounds = { 30*[1;1],5*[ 8; 10], 5*[ 12;  10]};

for i=1:numModes
    m(i).X = Polyhedron([eye(m(i).nx()); -eye(m(i).nx())],...
        [stateUpperBounds{i}; -stateLowerBounds{i}]);
    m(i).U = Polyhedron([eye(m(i).nu()); -eye(m(i).nu())],...
        [inputUpperBounds{i}; -inputLowerBounds{i}]);
end

%============================= COST MATRICES ==============================
m(1).Q = 4*[3.6 -3.8; -3.8 4.87];
m(2).Q = [10  -3  ; -3   8   ];
m(3).Q = [ 5  -4.5; -4.5 5   ];
m(3).Q = [ 10  -4.5; -4.5 5   ];

m(1).R = 2.6;
m(2).R = 1.165;
m(3).R = 1.111;

%LQR Calculation on each mode
arrayfun(@(x) x.idare(), m);