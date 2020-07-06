function m = maxOverSet(fun, Set)
m = -Inf;
if(nargin < 2)
    t = 0:pi/50:2*pi;
    x = [cos(t); sin(t)];
    for i=1:size(x,1)
        m = max(fun(x(:,i)), m);
    end
else
    x = Set.V;
    for i=1:size(x,1)
        m = max(fun(x(i,:)'), m);
    end
end
end