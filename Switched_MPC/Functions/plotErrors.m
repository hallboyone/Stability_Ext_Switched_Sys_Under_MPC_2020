function plotErrors(X, Y, costs)
figure
for i=1:size(costs, 1)
    names = fieldnames(costs{1});

    for name = 1:size(names, 1)
        if all(names{name} =="low") ||...
                all(names{name} =="true") ||...
                all(names{name} =="norm") ||all(all(isnan(costs{i}.(names{name}))))
            continue;
        end
        subplot(size(costs, 1), 1, i);
        
        err = (costs{i}.(names{name}) - costs{i}.true)./costs{i}.true;
        mx = max(max(err));
        surf(X, Y, err, 'FaceColor', 'interp')
        caxis([-.1, mx])
        zlim([-0.1, max(2*mx, 1)])
        title("Error using " + names{name},'Interpreter','none')
        
        [err_max, I1] = max(err);
        [err_max, I2] = max(err_max);
        I1 = I1(I2);
        text(X(I2), Y(I1), err_max+.25, "  Max error: "+ num2str(err_max))
    end
end
end