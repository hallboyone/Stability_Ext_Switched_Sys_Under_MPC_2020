function sandbox(c, m, ops)
% %% SYSTEM SETUP
% gamma = zeros(1,numel(mode)); %s.t. J^*(t + 1) < gamma * J^*(t)
% 
% %% Compute Gamma
% fprintf("\nTesting tighter gamma bounds\n");
% tic
% 
% for m_i = 1:numel(mode)
%     %======= Compute the smallest step down while dwelling in this mode ========
%     [costSet, fragments, Ksets] = fragmentQuadCost(mode(m_i), mode(m_i).S(end));
%     for frag_i = 1:numel(fragments)
%         gamma(m_i) = max(gamma(m_i), 1+frankWolfe(-(mode(m_i).Q+Ksets{frag_i}'*mode(m_i).R*Ksets{frag_i}), costSet{frag_i}, fragments(frag_i).V));
%     end
% end
% toc
% 
% X_vals = ops.mesh.value(1:ops.mesh.numPoints);
% 
% gamma_sim = zeros(1,numel(mode));
% 
% for m_i = 1:numel(mode)
%     gamma_grid = nan(ops.mesh.sz);
%     u = c.Inputs{m_i};
%     for x_i=1:ops.mesh.numPoints
%         if(~isnan(u(x_i)))
%             x = X_vals(:,x_i);
%             if(norm(x, 2) > 0.25)
%                 x_plus = mode(m_i).f(x, u(x_i));
%                 gamma_grid(x_i) = mode(m_i).MPCCost(x_plus)/mode(m_i).MPCCost(x);
%                 gamma_sim(m_i) = max(gamma_sim(m_i), gamma_grid(x_i));
%             end
%         end
%     end
%     figure
%     surf(ops.mesh.vals(1,:), ops.mesh.vals(2,:), gamma_grid)
%     disp(gamma_sim(m_i)/gamma(m_i))
% end
mdt = [15;21;13];
costs = zeros(20,30);
figure;
hold on
for sample = 1:20
    x = [7; 6.5];
    st = 0;
    m_i = randi(3);
    for step = 1:30
        costs(sample, step) = m(1).MPCCost(x);
        x = m(m_i).fMPC(x);
        if st==mdt(m_i)
            new_m_i = randi(3);
            while new_m_i==m_i
               new_m_i = randi(3);
            end
            m_i = new_m_i;
            st = 0;
        end
    end
    plot(1:30, costs(sample,:))
end

end
