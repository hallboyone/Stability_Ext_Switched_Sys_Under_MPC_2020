function brkt_w = brkt(w)
%BRKT  Flips between the vector and skew symmetric representations
%   If the value of w is a vector, then this function generates the
%   skew symmetric form of either a 1x3 axis of a 1x6 twist. If it is a
%   skew symmetric matrix, then it computes the coorisponding vector. A
%   skew symmetric matrix must satisfy A=-A'.

%Figure out what dir we are converting
dim = size(w);

%Converting to an axis
if all(dim == [3,3])
    if any(w ~= -(w'))
        fprintf("Not skew symmetric\n")
        return;
    end
    brkt_w = [w(3,2); w(1,3); w(2,1)];
    
%Converting to a twist
elseif all(dim == [4,4])
    if any(any(w(1:3,1:3) ~= -(w(1:3,1:3)')))||...
            any(w(4, :)~=0)
        fprintf("Not valid\n")
        return;
    end
    brkt_w = [brkt(w(1:3, 1:3)); w(1:3,4)];
    
%Converting from an axis
elseif all(dim == [3,1])||all(dim == [1,3])
    brkt_w = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
    
%Converting from a twist
elseif all(dim == [6,1])||all(dim == [1,6])
    %make sure the vector is upright
    if dim(2)>dim(1)
        w = w';
    end
    brkt_w = [brkt(w(1:3)), w(4:6); 0 0 0 0];
else
    fprintf("Invalid input\n")
    brkt_w = [];
end

return
end