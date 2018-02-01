% Author:   Nicole Chan
% Created:  1/11/18
% Description: Implementation for a cell-array-based list. 
%
% NOTE: cell-arrays: (i) references the i-th cell element, whereas {i}
% dereferences the object wrapped inside the cell
%
classdef listObj < handle
    properties
        items = listObj.empty;  % objects in list
        last  = [];             % index of last object in list (equivalent to size of list)
    end
    
    methods
        function lst = listObj(obj)
            if nargin == 1
                if iscell(obj)
                    dim = numel(obj);
                    lst.items = obj;
                    lst.last = dim;
                else
                    lst.items = {obj};
                    lst.last = 1;
                end
            else
                lst.items = {};
                lst.last = 0;
            end
        end
        
        function push(lst,obj)
            if iscell(obj)
              dim = numel(obj);
              lst.items(lst.last+1:lst.last+dim) = obj;
              lst.last = lst.last+dim;
            else
              lst.last = lst.last+1;
              lst.items(lst.last) = {obj};
            end
        end
        
        % Remove the element at index ind; if ind not specified, removes
        % last element
        function pop(lst,ind)
            if nargin==1
                if lst.last~=0
                    lst.items(lst.last) = []; % use 'lst.items{lst.last} = []' if not deleting cell element
                    lst.last = lst.last-1;
                else
                    disp('Empty list, nothing to pop.');
                end
            elseif nargin==2 && ind > 0 && ind <= lst.last
                lst.items(ind) = [];
                lst.last = lst.last-1;
            else
                error('Invalid input arguments.');
            end
        end
        
        function y = isempty(lst)
            y = (lst.last == 0);
        end
        
        function y = size(lst)
            y = lst.last;
        end
        
        function y = getLast(lst)
            if lst.last~=0
              y = lst.items{lst.last}; 
            else
              y = [];
            end
        end
        
        function y = getItem(lst,ind)
            if nargin==2 && ind > 0 && ind <= lst.last
                y = lst.items{ind};
            else
                error('Invalid input arguments.');
            end
        end
    end % end methods
end % end classdef