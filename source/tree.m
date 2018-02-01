% Author:   Nicole Chan
% Created:  1/10/18
% Description: The 'tree' object is 
%
classdef tree < handle
    properties
        rootNode@treeNode   % handle to root node
        size                % number of nodes in tree
        height              % height of tree
        currNode@treeNode   % handle to some node
    end
    
    methods
        % Constructor
        function tree = tree(root)
            if nargin ~= 1
                tree.rootNode = treeNode();
                tree.size = 0;
                tree.currNode = [];
            else
                tree.rootNode = treeNode(root);
                tree.size = 1;
                tree.currNode = tree.rootNode;
            end
            tree.height = 0;
        end
        
        % Returns handle to root node
        function root = getRoot(treeObj)
            root = treeObj.rootNode;
        end
 
    end % end methods
end % end classdef