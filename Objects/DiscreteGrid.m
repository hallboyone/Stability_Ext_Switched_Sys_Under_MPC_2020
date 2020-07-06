% DISCRETEGRID Object to handle the indexing of a discrete grid.
% This object is used to sample points over a gridded space. Up to 4D is
% supported. 
%
% Construction : 
% mesh = DISCRETEGRID(edgePoints, dim) creates a mesh with the given
% diminsion with edges given by the edgePoints vector. Best for systems
% where the scale of each diminsion is relativly uniform.
%
% mesh = DISCRETEGRID(edgePointsCell) creates a mesh grid with a diminsion
% equal to the number of cells in edgePointsCell. Each cell contains the
% edge data for that diminsion in a vector. Best for systems where the
% scale of the dimisnions are very differant.
%
% Examples :
% mesh = DISCRETEGRID(-10:.1:10, 3) Create a 3D mesh grid with 8120601
% points
%
% mesh = DISCRETEGRID({-10:.1:10, -1:.01:1}) Create a 2D mesh rangeing from
% -10 to 10 in the first diminsion and -1 to 1 in the second diminsion.

classdef DiscreteGrid
    properties
        vals; %Matrix or vector of edge values
        dim;  %Diminsion of grid 
        sz;   %Vector holding the number of points in each diminsion
        numPoints; %Total number of points in the grid
    end
    
    methods
        function obj = DiscreteGrid(vals, dim)
            %Construct a DiscreteGrid object
            
            %Make sure vals is 1D array or cell
            if min(size(vals)) > 1
                error("vals should be 1D array or cell array");
            end
            
            %If diminsion was provided. Use single vals vector
            if exist('dim', 'var')
                %Check the diminsion
                if dim > 4
                    error("The discreteGrid object only works up to 4 diminsions");
                end
                
                %If the user passed in a single element cell array, convert
                %to vector
                if iscell(vals)
                    if numel(vals)==1 && min(size(vals{1}))==1
                        vals = cell2mat(vals);
                    else
                        error("If specifying diminsion, vals should be a 1D array");
                    end
                end
                
                %Make a sz array 
                obj.sz = repmat(max(size(vals)), 1, dim);
                
                %Compute the number of points
                obj.numPoints = obj.sz(1)^dim;
                
                %Make sure the grid array won't be too large
                if ispc
                    obj.checkMem();
                end
                
                %Same the values and dim
                obj.vals = repmat(vals, dim, 1);
                obj.dim = dim;

            else%, use cell array of edgeVals
                
                %Check diminsion
                if numel(vals) > 4
                    error("The discreteGrid object only works up to 4 diminsions");
                end
                
                %If the vals are not a cell array, make sure they are a 1D
                %vector for a 1D system
                if ~iscell(vals)
                    if min(size(vals))==1
                        vals = {vals};
                    else
                        error("Vals must be a cell array of 1D vectors");
                    end
                end
                
                %Save the diminsion
                obj.dim = numel(vals);
                
                %Iterate through each element of the cell array and get their
                %size, save their values, and update the number of points.
                obj.sz = zeros(1,obj.dim);
                obj.numPoints = 1;
                
                for i=obj.dim:-1:1
                    obj.sz(i) = max(size(vals{i}));
                    obj.numPoints = obj.numPoints * obj.sz(i);
                end
                
                obj.vals = nan(obj.dim, max(obj.sz));
                for i=obj.dim:-1:1
                    obj.vals(i, 1:obj.sz(i)) = vals{i};
                end
                
                %Make sure the DiscreteGrid does not take up to much space
                if ispc
                    obj.checkMem();
                end
            end
            
            %Set things very close to 0 exactly to 0
            obj.vals(abs(obj.vals)<10e-15) = 0;
        end
            
        function val = value(obj, i)
            %Returns the value of the grid at linear index i
            if obj.dim==1
                val = obj.vals(i);
            elseif obj.dim==2
                [ix, iy] = ind2sub(flip(obj.sz), i');
                val = [obj.vals(1, iy); obj.vals(2, ix)];
            elseif obj.dim==3
                [ix, iy, iz] = ind2sub(flip(obj.sz), i');
                val = [obj.vals(1, iz); obj.vals(2, iy); obj.vals(3, ix)];
            else
                [i1, i2, i3, i4] = ind2sub(flip(obj.sz), i');
                val = [obj.vals(1, i4); obj.vals(2, i3); obj.vals(3, i2); obj.vals(3, i1)];
            end     
        end
        
       
        function areEqual = equals(obj, rhs)
            %Checks if two DiscreteGrid objects are equal
            areEqual = (obj.dim == rhs.dim) && (isequal(obj.vals(~isnan(obj.vals)), rhs.vals(~isnan(rhs.vals))));
        end
         
        function idx = origin(obj)
            %Returns the linear index of the origin
            
            %Get the indicies of the columns containing 0
            [~, zero_idx] = find(obj.vals==0 & ~isnan(obj.vals));
            
            %If at least one column does not contain 0, mesh does not contain 0
            if numel(zero_idx) ~= obj.dim
                warning("Grid does not contain the origin");
                idx = [];
            else
                %Convert the sub indexes computed above to a linear one.
                if obj.dim==1
                    idx = sub2ind(obj.sz, zero_idx);
                elseif obj.dim==2
                    idx = sub2ind(flip(obj.sz), zero_idx(2), zero_idx(1));
                elseif obj.dim==3
                    idx = sub2ind(flip(obj.sz), zero_idx(3), zero_idx(2), zero_idx(1));
                else
                    idx = sub2ind(flip(obj.sz), zero_idx(4), zero_idx(3), zero_idx(2), zero_idx(1));
                end
            end
        end
        
        function surf(obj, Z, plotMesh)
            %Plots the surface Z over a 2D grid object
            if obj.dim==2
                surf(obj.vals(1,1:obj.sz(1)), obj.vals(2,1:obj.sz(2)), reshape(Z, flip(obj.sz)));
                if exist('plotMesh', 'var')
                    if plotMesh
                        hold on
                        mesh(obj.vals(1,1:obj.sz(1)), obj.vals(2,1:obj.sz(2)), zeros(flip(obj.sz)));
                    end
                end
            else
                warning("Can only make a surf plot of 2D grids")
            end
        end
    end
    
    methods (Access = private)

        function checkMem(obj)
            %Make sure that we will not have more points than computer can
            %handle
            mem = memory;
            largestArray = mem.MaxPossibleArrayBytes/8;
            if obj.numPoints > largestArray
                error("Grid has more points than max array size");
            elseif obj.numPoints > largestArray/4
                warning('Grid has a large number of points compared with the largest possible arrray. %d vs %d. Performance may suffer',...
                    largestArray, obj.numPoints);
            end
        end
    end
end