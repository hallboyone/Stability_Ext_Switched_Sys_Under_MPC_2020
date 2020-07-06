function minDwellTimes(sCost, gamma, name)
% MINDWELLTIMES  Computes the min dwell times given the cost of a switch
% and the decrease at each time step.

D = ones(size(sCost));
for i = 1:numel(gamma)
    for j = 1:numel(gamma)
        if i ~= j
            D(i,j) = max(1, ceil(-log(sCost(i, j)) / log(gamma(i))));
        end
    end
end

fprintf("\nThe MM-MDTs using the " + name + " method are:\n");
for i = 1:numel(gamma)
    fprintf(" - ")
    for j = 1:numel(gamma)
        if i ~= j
            fprintf(" - D(%d,%d) = %6d    ", i, j, D(i,j));
        end
    end
    fprintf("\n")
end

fprintf("\nThe M-MDTs using the " + name + " method are:\n");
%For each mode
for i=1:numel(gamma)
    disp(" - " + max(D(i,:)) + " for Mode " + i);
end

fprintf("\nThe MDT using the " + name + " method is %d\n", + max(max(D)));
end