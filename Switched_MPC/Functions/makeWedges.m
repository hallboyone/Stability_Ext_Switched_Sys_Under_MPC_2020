function [wedges, vertexIndex] = makeWedges(S)

%Ensure the polyhedron is in its min rep
S.minHRep;
S.minVRep;

%Save the set of inequalities and verticies
A = S.A;
b = S.b;
V = S.V;
dim = S.Dim;

%Row i is inequality i, col j is vertex j. Element (i,j) is true if vertex j
%activivates inequality i. False else.
adjMat = false(numrow(A), numrow(V));
for i=1:numrow(A)
    adjMat(i,abs(A(i,:) * V' - b(i)) < 10e-9) = true;
end

%The i'th row contains the vertex indicies of inequality i. They will be sorted
%in the next for loop so that each adjacent index in the matrix shares an edge.
Vidx = nan(numrow(A), numrow(V));

%Iterate through each inequality
for H_i=1:numrow(A)
    %Make a flag of the facets which could be neighbors (all but the current)
    potentialNeighbors = true(numrow(A), 1);
    potentialNeighbors(H_i) = false;
    
    %Make an unordered list of verticies in facet. These are what will be sorted
    %and added to the H_i row of Vidx
    unorderedV = find(adjMat(H_i,:));
    
    %Save the first vertex as the starting point in the list of sorted verticies
    %and remove it from the list of unsorted
    Vidx(H_i,1) = unorderedV(1);
    unorderedV(1) = [];
    
    %Variable of the number unsorted vertices in facet
    numV = numel(unorderedV);
    
    %Make a chain a verticies around the facet
    for V_i = 1:numV-1
        %Facets that contain the latest vertex added to the sorted array
        sharedHIdx = find(adjMat(:, Vidx(H_i,V_i)) & potentialNeighbors)';
        
        %For each neighbor facet, check if if contains an element in the
        %unorderedV. If so, that is the neighbor point
        for k=sharedHIdx
            %Used this neighbor, no longer an option
            potentialNeighbors(k) = false;
            
            neighbor_contains = adjMat(k, unorderedV);
            %Does the suspect facet contain any other points in current
            %facet
            if any(neighbor_contains)
                %Add the sorted vertex
                Vidx(H_i,V_i+1) = unorderedV(neighbor_contains);
                
                %Remove the sorted vertex from the list of unsorted
                unorderedV(neighbor_contains) = [];
                break;
            end
        end
    end
    %If the for loop was never run
    if isempty(V_i)
        V_i = 0;
    end
    Vidx(H_i,V_i+2) = unorderedV;
end

%Get the total number of wedges and make vars to hold them and their vertex
%index
numWedges = sum(sum(~isnan(Vidx), 2)-dim+1);
vertexIndex = nan(numWedges, dim);
wedges = repmat(Polyhedron, numWedges, 1);

% %Triangulate each facet: Greedy algorithm
% figure
% w_i = 1;
% for i = 1:numrow(A)
%     V_face = V(Vidx(i,~isnan(Vidx(i,:))),:);
%     idx = optimalWedges(V_face, 1, numrow(V_face));
%     for k = 1:numrow(idx)
%         wedges(w_i) = Polyhedron([V_face(idx(k,:), :); zeros(1,dim)]);
%         w_i = w_i + 1;
%     end
% end
% plot(wedges);

%Triangulate each facet: Fan algorithm
H_i = 1;
V_i = 2;
for i = 1:numWedges
    %Move to the next row if needed
    if isnan(Vidx(H_i, V_i+dim-2))
        H_i = H_i+1;
        V_i = 2;
    end
    
    vertexIndex(i, 1) = Vidx(H_i, 1);
    for col=1:(dim-1)
        vertexIndex(i, col+1) = Vidx(H_i, V_i+col-1);
    end
    wedges(i) = Polyhedron([V(vertexIndex(i, :), :); zeros(1,dim)]);
    V_i = V_i + 1;
end
% plot(wedges);
end

function idx = optimalWedges(V, i, j)
if j<=i+1
    idx  = [];
else
    cost = inf;
    for k = (i+1):(j-1)
        c3 = sqrt((V(i,:) - V(j,:))*(V(i,:) - V(j,:))') +...
             sqrt((V(i,:) - V(k,:))*(V(i,:) - V(k,:))') +...
             sqrt((V(j,:) - V(k,:))*(V(j,:) - V(k,:))');
        if (cost > c3)
            opt_k = k;
        end
    end
    idx = [optimalWedges(V, i, opt_k); optimalWedges(V, opt_k, j); i, opt_k, j];
end
end
