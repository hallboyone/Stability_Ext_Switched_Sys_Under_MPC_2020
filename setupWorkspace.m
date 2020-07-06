function setupWorkspace()
%Add all the files in the current folder to path
disp("Adding folders to path...")

%Array of all the string we wish to ignore
patterns2Ignore = '.git';

%Make a cell array of all possible paths to folders
allPaths = split(genpath('.'), [";",":"]);

%Remove empty final cell
allPaths(end) = [];

%Iterate through possible paths. Add them if they do not contain any patterns to
%ignore
for pathIdx = 1:numel(allPaths)
    path = allPaths{pathIdx};
    if ~contains(path, patterns2Ignore)
        addpath(path);
    end
end

disp("Setup complete!");
end