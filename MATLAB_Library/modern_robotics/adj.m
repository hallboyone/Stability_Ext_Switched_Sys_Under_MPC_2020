function adj_T = adj(T)
%ADJ  Generate the adjoint form of a transform matrix
%   Function to compute the 6x6 adjoint matrix of a given transform matrix
%   T. The adjoint form is used to change the reference axis of a twist.

R = T(1:3, 1:3);
p = T(1:3, 4);
adj_T = [R, zeros(3, 3); brkt(p)*R R];
end