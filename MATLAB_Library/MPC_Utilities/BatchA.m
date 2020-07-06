function A_batch = BatchA(A, N)
% BATCHA - Returns the A matrix used to describe a series of LTI system
% updates as a single time step.

A_batch = zeros(size(A,1) * (N+1), size(A, 2));

for i=0:N
    A_batch((i*size(A,1)+1):(i+1)*size(A,1),:) = A^i;
end

end