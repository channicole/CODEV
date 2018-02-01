% Author:   Nicole Chan
% Created:  1/10/18
% Description: The 'treeNode' object is 
%
classdef treeNode < handle
    properties
        data            % main node info
        parent          % handle to parent node
        children        % array of children node handles (for now, using custom implementation of list)
    end
    
    methods
        % Constructor
        function node = treeNode(data,parent,children)
            if nargin < 3 || ~isa(children,'listObj')
                children = listObj();
            end
            if nargin < 2
                parent = [];
            end
            if nargin < 1
                data = [];
            end
            node.data = data;
            node.parent = parent;
            node.children = children;
        end
        
        % Check if this node is the root
        function flag = isRoot(node)
            flag = 0;
            if isempty(node.parent)
                flag = 1;
            end
        end
        
        % Check if this node is a leaf
        function flag = isLeaf(node)
            flag = 0;
            if node.children.size==0
                flag = 1;
            end
        end
                
        % Add child
        function node = insertChild(node,child)
            if isa(child,'treeNode')
                child.parent = node;
                node.children.push(child);
            else
                error('Input argument not of type treeNode.');
            end
        end
        
        % Remove child
        function node = removeChild(node,ind)
            if nargin==1
                node.children.pop();
            elseif nargin==2
                node.children.pop(ind); 
            end
        end
    end % end methods
end % end classdef