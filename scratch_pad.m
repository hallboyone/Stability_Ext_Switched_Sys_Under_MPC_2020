% c1 = [69,26,74]/255;
% c2 = [238,68,47]/255;
% c3 = [99,172,190]/255;
% c4 = [249,244,236]/255;
% p = panel();
% p.pack(1,2);
% 
% p(1,1).select()
% plot(m(1).T, 'color', c1)
% hold on
% plot(m(2).T, 'color', c2)
% plot(m(3).T, 'color', c3)
% plot(m(1).S(end), 'alpha', 0, 'linestyle', '--', 'linewidth', 1.5, 'edgecolor', c1)
% plot(m(2).S(end), 'alpha', 0, 'linestyle', '--', 'linewidth', 1.5, 'edgecolor', c2)
% plot(m(3).S(end), 'alpha', 0, 'linestyle', '--', 'linewidth', 1.5, 'edgecolor', c3)
% 
% legend(['Mode 1 - Terminal Set'; 
%         'Mode 2 - Terminal Set'; 
%         'Mode 3 - Terminal Set';
%         'Mode 1 - Feasable Set';
%         'Mode 2 - Feasable Set';
%         'Mode 3 - Feasable Set';],...
%         'Location', 'SE')
% title("Feasible and Terminal Sets");
% 
% p(1,2).select()
% [W, ~] = makeWedges(m(1).S(end));
% plot(W, 'edgecolor', c1, 'color', c4);
% title("Fragmented feasible set of mode 1")
% p.de.margin = 7;
% p.margin = [7 7 2 7];
% set(gcf,'position',[452   662   788   420])
% saveFig('fig_1', 'high')
% clf

p = 0.05;
a = 0;
steps = 500;
a_vals = nan(steps,1);
p_vals = nan(steps,1);
for n=1:steps
    a_vals(n) = a;
    p_vals(n) = p*(1-p)^(n-1);
    a = a + n*p*(1-p)^(n-1);
end
plot(a_vals)
figure
plot(p_vals)