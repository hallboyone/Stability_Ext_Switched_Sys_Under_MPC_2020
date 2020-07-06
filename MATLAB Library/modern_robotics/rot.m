function R = rot(axis, theta)
%ROT  Generate rotation matrix
%   Function to compute the 3x3 rotation matrix given a axis of rotation and an
%   angle (in radians) to rotate. The axis can either be given as a 1x3
%   vector, or a string ('x', 'y', or 'z') specifiying a coordinate axis
%
%   R = ROT('x', pi/2) generates a 90 deg rotation matrix about the x axis.
%
%   R = ROT([1, 2, 3], 1.23) generates a 1.23 rad rotation matrix about the given vector

if ischar(axis) %If the axis is specified as 'x', 'y', or 'z'
    switch axis
        case 'x'
            axis = [1 0 0];
        case 'y'
            axis = [0 1 0];
        case 'z'
            axis = [0 0 1];
        otherwise
            fprintf("Not a valid axis. Enter a 3x1 vector or either 'x', 'y', or 'z'\n")
            return
    end
end

%Norm the axis
axis = axis/norm(axis);

R = zeros(3);

R(1,1) = cos(theta)+axis(1)^2*(1-cos(theta));
R(2,2) = cos(theta)+axis(2)^2*(1-cos(theta));
R(3,3) = cos(theta)+axis(3)^2*(1-cos(theta));

R(1,2) = axis(1)*axis(2) * (1-cos(theta)) - axis(3)*sin(theta);
R(1,3) = axis(1)*axis(3) * (1-cos(theta)) + axis(2)*sin(theta);

R(2,1) = axis(1)*axis(2) * (1-cos(theta)) + axis(3)*sin(theta);
R(2,3) = axis(2)*axis(3) * (1-cos(theta)) - axis(1)*sin(theta);

R(3,1) = axis(1)*axis(3) * (1-cos(theta)) - axis(2)*sin(theta);
R(3,2) = axis(2)*axis(3) * (1-cos(theta)) + axis(1)*sin(theta);

end