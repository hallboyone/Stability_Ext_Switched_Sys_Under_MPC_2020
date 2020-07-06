function saveFig(name, quality, fig)
% SAVEFIG : Save the current figure or the figure passed as fig to the file
% ./figures/[name].png. The quility can be set to [1||'low'], [2||'med'],
% [3||'high'], or [4||'max']. Default is [2||'med']
% EXAMPLE 1: saveFig('my_fig', 'high', fig1);
% EXAMPLE 2: saveFig('my_other_fig');


if(nargin == 3)
    figure(fig);
end
name = convertStringsToChars(name);
if(nargin >= 2)
    switch quality
        case {1, 'low'}
            res = '-r75';
        case {2, 'med'}
            res = '-r150';
        case {3, 'high'}
            res = '-r300';
        case {4, 'max'}
            res = '-r450';
        otherwise
            res = '-r150';
    end
else
    res = '-r150';
end

if ~exist('./figures', 'dir')
    mkdir('./figures')
end

print(['./figures/', name], '-dpng', res);
end
        
        