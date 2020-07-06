function [batch_Q, batch_R] = BatchCost(Q, R, P, N)
batch_Q = blkdiag(kron(eye(N), Q), P);
batch_R = kron(eye(N), R);
end