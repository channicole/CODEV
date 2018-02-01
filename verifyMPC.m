% Author:   Nicole Chan
% Created:  1/23/18
% Description: Implements C2E2's algorithm for a specific class of hybrid
% models, which use the (linear) MPC formulated using a mp-LP solver. 
% The model transitions periodically (as determined by sampling size in
% MPC) and the number of discrete modes correspond to the number of regions
% in the MPC solution.
% The controller is linear state feedback, but the plant dynamics may be
% nonlinear.
% For now, we only "refine covers" of the initial set, not subsequent
% reachsets. Furthermore, we do not search for counterexamples and decouple
% this function from safety certification.
%
function [reach,safe,Thorizon]=verifyMPC(inputHandle)
% dbstop if error
addpath(genpath('source'));
%% Setup
safe = 0;                           % 0: UNKNOWN, 1: SAFE

% Load user input from the file in the input arg
if ~isa(inputHandle,'function_handle')
    error('Input should be a function handle to the input file.');
end
inParams = inputHandle();           % instantiates function handles to each input entity needed from user
inFunNames = fieldnames(inParams);  % get function handle names from user input, stored in a cell array
if length(inFunNames) < 6           % TODO: update with inputParams
    error('Incorrect number of function handles. Follow template for the inputParams() method.');
else
    template = {'flowEq';'modelJac';'unsafeStates';'initStates';'vPar';'MPCprob'}; % TODO: update with inputParms
    if ~isempty(setdiff(template,inFunNames))
        error('Incorrect function handle names. Follow template for the inputParams() method.');
    end
end

%% Compute the MPC solution if one isn't already specified
if ~isempty(setdiff('MPCsol',inFunNames))
    [sol.Pn,sol.Fi,sol.Gi,~,~,~] = MPCsolve(inParams.MPCprob);
    inParams.MPCsol = sol; % add the MPC solution to the inParams struct
end

%% Instantiate unsafeSets if it is dependent on MPC solution
if isa(inParams.unsafeStates,'function_handle')
    inParams.unsafeStates = inParams.unsafeStates(length(inParams.MPCsol.Pn));
end

%% Safety Verification: Generate covers using DFS
% Get the time horizon desired for the reachtube
Thorizon = (inParams.vPar.Thorizon(1):inParams.vPar.simStep:inParams.vPar.Thorizon(end))';

% Check if the cover intersects with multiple discrete mode invariants
[x0,rad0,ball0,Theta0,MPC0] = poly2ball(inParams.initStates,inParams.MPCsol.Pn);
numRegions = size(rad0,1);

if numRegions > 1
    % Initialize tree of coverObj-objects
    covTheta = coverObj(Theta0,[],[],[],[],Thorizon,[]);
    CTree = tree(covTheta);
    
    % Add children covers with single mode covers
    for i=1:numRegions
        % Get MPC solution for each region
        MPC.Fi = inParams.MPCsol.Fi{MPC0(i)};
        MPC.Gi = inParams.MPCsol.Gi{MPC0(i)};
        MPC.ind = MPC0(i);
        % Create coverObj and add to tree
        CTree.rootNode.insertChild(treeNode(coverObj(Theta0(i),x0(i,:)',...
            x0(i,:),MPC,ball0(i),Thorizon,rad0(i,:)'),CTree.rootNode));
    end
    CTree.height = CTree.height + 1;
    CTree.size = CTree.size + numRegions;
    CTree.currNode = CTree.currNode.children.getItem(1);
elseif numRegions == 1 
    % Get MPC solution for each region
    MPC.Fi = inParams.MPCsol.Fi{MPC0};
    MPC.Gi = inParams.MPCsol.Gi{MPC0};
    MPC.ind = MPC0;
    % Initialize tree of coverObj-objects
    covTheta = coverObj(Theta0,x0,x0',MPC,ball0,Thorizon,rad0);
    CTree = tree(covTheta);
else
    error('Check MPC solution and initial set of states.');
end

% Compute the reachtubes for each corresponding 'cover' in the stack
while ~isempty(CTree.currNode)
    switch CTree.currNode.data.isSafe
        % Quit program
        case -2
            CTree.currNode = treeNode.empty();
            disp('Maximum number of partitions reached. Safety result unknown.');
            
        % Compute reachtube for the cover and check its safety
        case -1
            CTree.currNode = computeReach(CTree.currNode,inParams);
            
        % Partition the cover
        case 0
            if CTree.height == inParams.vPar.maxPart
                CTree.currNode.data.setSafeFlag(safetyEnum.ExitPartitionBnd);
            else
                newCov = CTree.currNode.data.partitionCover();
                for i=1:length(newCov)
                    CTree.currNode.insertChild(treeNode(newCov{i},CTree.currNode));
                end
                % Update currNode pointer to left-most child
                CTree.currNode = CTree.currNode.children.getItem(1);
                % Update tree parameters
                CTree.height = CTree.height + 1;
                CTree.size = CTree.size + length(newCov);
            end
            
        % Update currNode to next cover that needs to be checked
        case 1
            % Check if root
            if isempty(CTree.currNode.parent)
                reach = CTree.currNode.data;
                CTree.currNode = treeNode.empty();
                break;
            end
            
            % Check for left sibling with safeFlag==1 
            lSib = CTree.currNode.parent.children.getItem(1);
            if lSib~=CTree.currNode
                if lSib.data.isSafe~=1
                    error('Left sibling cover should have had reachtube computed to be RobustlySafe.');
                end
                % Combine covers that have been computed in currNode and left sibling
                CTree.currNode.data.coverUnion(lSib.data);
                % Pop left sibling
                CTree.currNode.parent.removeChild(1);
                CTree.size = CTree.size - 1;
            end
            
            % Check for existence of a right sibling
            if CTree.currNode.parent.children.size > 1
                CTree.currNode = CTree.currNode.parent.children.getItem(2);
            else
                % If no siblings left to check, then update parent:
                % Copy child level to parent
                CTree.currNode.parent = CTree.currNode;
                % Update currNode to parent 
                CTree.currNode = CTree.currNode.parent;
                % Delete child
                CTree.currNode.removeChild(1);
                CTree.size = CTree.size - 1;
                CTree.height = CTree.height - 1;
                % Update safety flag
                CTree.currNode.setSafeFlag('RobustSafe');
            end
            
        otherwise
            error('Not sure how you got here. This parameter is enumerated and should automatically be checked/generate an error.');
    end
end

% Return reachtube
[~,~,reach] = CTree.rootNode.data.reachtube.getProperties();
safe = 1;
end