function idxSet = getIndexSet(dim, n)
%GETINDEXSET Generates a (d x n^d) array of the grid indicies in d-D
%   Ceates the array of indices, starting with [1,1,...,1] in the first
%   row, [2,1,...,1] in the second, etc, to indicate points in a gridded
%   dim diminsional space. Since the output can be very large, the output
%   is actually an anonymous function which generates the row idx when
%   called ie, idxSet(idx)
%EX - getIndexSet(3, 2)
% Returns [1, 1, 1;
%          2, 1, 1;
%          1, 2, 1;
%          2, 2, 1;
%          1, 1, 2;
%          2, 1, 2;
%          1, 2, 2;
%          2, 2, 2];

%idx = mod(floor([(i-1)./(n.^[0:(d-1)])]),n) + 1;

idxSet = [1:n]';
idx_cur = idxSet;
for i=2:dim
    idx_cur = sort(repmat(idx_cur, n, 1));
    idxSet = [idx_cur, repmat(idxSet, n, 1)];
end
flip(idxSet, 2);
end
    
