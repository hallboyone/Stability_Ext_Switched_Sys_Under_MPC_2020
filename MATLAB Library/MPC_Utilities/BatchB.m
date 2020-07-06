function B_batch = BatchB(A, B, N)
% BATCHB - Returns the B matrix used to describe a series of LTI system
% updates as a single time step.

B_batch = zeros(size(B,1)*(N+1), size(B,2)*N);

for i=1:N+1
    for j=1:i-1
       B_batch((i-1)*size(B,1)+1:i*size(B,1), (j-1)*size(B,2)+1:j*size(B,2)) = A^(i-j-1) * B;
    end
end

end