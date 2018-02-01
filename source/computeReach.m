% Author:   Nicole Chan
% Created:  1/11/18
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a hybrid system.
%
% When this function is called, the input 'cover' should already correspond
% to the dynamics of no more than one discrete mode.
%
function cover=computeReach(cover,in)
    if nargin~=2
        error('Not enough input arguments.');
    end
    %% Start building tree of reachObj nodes (these nodes are eventually concatenated into coverObj.Reach)
    RTree = tree(cover.data.getReachtube); % Note RTree.currNode is set to the root node
    
    while ~isempty(RTree.currNode)
        % Compute Post (only leaf nodes do NOT need Post)
        if RTree.currNode.data.needsPost()==1 % RTree.currNode.data should be type 'reachtubeObj'
            % If transitions are only checked periodically every deltaStep
            if in.vPar.deltaStep > in.vPar.simStep
                tf = min(RTree.currNode.data.T(end),RTree.currNode.data.T(1)+in.vPar.deltaStep);
                % If the unsafe set is dependent on the mode/MPC sol
                if isfield(in.unsafeStates, 'funHand')
                    % Check if the current unsafe set has been instantiated before
                    regInd = RTree.currNode.data.MPCi.ind;
                    if isempty(in.unsafeStates.array{regInd})
                        in.unsafeStates.array{regInd} = in.unsafeStates.funHand(...
                            RTree.currNode.data.MPCi.Fi,RTree.currNode.data.MPCi.Gi);
                    end
                    % Compute post
                    [outReach,outSafe] = computePost(RTree.currNode.data.Theta,...
                        RTree.currNode.data.T(1:find(abs(RTree.currNode.data.T-tf)<=eps,1)),...
                        in.flowEq,in.MPCsol,RTree.currNode.data.MPCi,in.modelJac,...
                        in.unsafeStates.array{regInd},in.vPar.epsilonConst,in.vPar.LipConst,1);
                % If the unsafe set is static
                else
                    [outReach,outSafe] = computePost(RTree.currNode.data.Theta,...
                        RTree.currNode.data.T(1:find(abs(RTree.currNode.data.T-tf)<=eps,1)),...
                        in.flowEq,in.MPCsol,RTree.currNode.data.MPCi,in.modelJac,...
                        in.unsafeStates,in.vPar.epsilonConst,in.vPar.LipConst,1);
                end
            % If transitions are checked after every simStep
            else
                % If the unsafe set is dependent on the mode/MPC sol
                if isfield(in.unsafeStates, 'funHand')
                    % Check if the current unsafe set has been instantiated before
                    regInd = RTree.currNode.data.MPCi.ind;
                    if isempty(in.unsafeStates.array{regInd})
                        in.unsafeStates.array{regInd} = in.unsafeStates.funHand(...
                            RTree.currNode.data.MPCi.Fi,RTree.currNode.data.MPCi.Gi);
                    end
                    % Compute post
                    [outReach,outSafe] = computePost(RTree.currNode.data.Theta,...
                        RTree.currNode.data.T,in.flowEq,in.MPCsol,RTree.currNode.data.MPCi,...
                        in.modelJac,in.unsafeStates.array{regInd},in.vPar.epsilonConst,...
                        in.vPar.LipConst,0);
                    
                % If the unsafe set is static
                else
                    [outReach,outSafe] = computePost(RTree.currNode.data.Theta,...
                        RTree.currNode.data.T,in.flowEq,in.MPCsol,RTree.currNode.data.MPCi,...
                        in.modelJac,in.unsafeStates,in.vPar.epsilonConst,in.vPar.LipConst,0);
                end
                % Check safety of the computed reachtube
                if outSafe==-1
                    % Get the part of the reachtube that intersects with the mode invariant
                    [~,~,newReach] = outReach.getProperties();
                    [newReach,nextTheta,nextT,nextMPCind] = getInvariantIntersect(newReach,in.MPCsol.Pn);
                    
                    % Update outReach object
                    outReach.updateReach(newReach);
                    outReach.T = outReach.T(1:length(newReach));
                    
                    outSafe = checkSafety(newReach,in.unsafeStates);
                end
            end
            
            % If found a possibly unsafe region, stop reachtube computation
            if ~outSafe
                cover.data.setSafeFlag(safetyEnum.NotSafe);
                % Update currNode
                RTree.currNode = treeNode.empty();
                return
                
            % If safe, continue computing reachtube
            elseif outReach.T(end) < RTree.currNode.data.T(end)
                % Update node with computed sub-reachtube
                [~,~,newReach] = outReach.getProperties();
                RTree.currNode.data.updateReach(newReach);
                
                % Add resulting covers to check for periodic guards
                if in.vPar.deltaStep > in.vPar.simStep
                    % Add nodes/sets of initial sets for next interval of time
                    newT = RTree.currNode.data.T(find(...
                        abs(RTree.currNode.data.T-outReach.T(end))<=eps,1):end);
                    try
                        [newx0,newrad,newBall,newTheta,MPCind] = poly2ball(newReach(end),in.MPCsol.Pn);
                    catch
                        % infeasible region reached
                        cover.data.setSafeFlag(safetyEnum.ExitInfeasible);
                        return
                    end
                    for i=1:size(newrad,1)
                        % Get MPC solution for each region
                        MPC.Fi = in.MPCsol.Fi{MPCind(i)};
                        MPC.Gi = in.MPCsol.Gi{MPCind(i)};
                        MPC.ind = MPCind(i);
                        % Create new reachtubeObj and add to tree
                        RTree.currNode.insertChild(treeNode(reachtubeObj(...
                            newTheta(i),newx0(i,:),[],MPC,newBall(i),newT,...
                            newrad(i,:)),RTree.currNode));
                        RTree.size = RTree.size+1;
                    end
                    
                % Add resulting covers to check for non-temporal guards
                else
                    for i=1:length(nextT)
                        % Get next time horizon
                        newT = RTree.currNode.data.T(find(...
                            abs(RTree.currNode.data.T-nextT(i))<=eps):end);
                        [newx0,newrad,newBall] = poly2ball(nextTheta(i));
                        % Get MPC solution for each region
                        MPC.Fi = in.MPCsol.Fi{nextMPCind(i)};
                        MPC.Gi = in.MPCsol.Gi{nextMPCind(i)};
                        MPC.ind = nextMPCind(i);
                        % Create new reachtubeObj and add to tree
                        RTree.currNode.insertChild(treeNode(reachtubeObj(...
                            nextTheta(i),newx0,[],MPC,newBall,newT,newrad)...
                            ,RTree.currNode));
                        RTree.size = RTree.size+1;
                    end
                end
                
                RTree.height = RTree.height+1;
                
                % Set currNode to point to the left-most child
                RTree.currNode = RTree.currNode.children.getItem(1); %RTree.currNode.children.items{1}
                
            % If safe and completed a branch of the reachtube computation
            else
                % Update node with computed sub-reachtube
                [~,~,newReach] = outReach.getProperties();
                RTree.currNode.data.updateReach(newReach);
                
                % If a left sibling exists
                lSib = RTree.currNode.parent.children.getItem(1);
                if lSib~=RTree.currNode
                    if lSib.data.needsPost()
                        error('Left sibling should have Post computed already.');
                    end
                    % Combine reachsets computed in these sibling nodes
                    RTree.currNode.data.reachUnion(lSib.data);
                    % Pop left sibling
                    RTree.currNode.parent.removeChild(1); % if lSib is referencing this node, then we might need to clear lSib or Matlab might not delete the object itself
                    RTree.size = RTree.size-1;
                elseif RTree.currNode.parent.children.size > 1
                    RTree.currNode = RTree.currNode.parent.children.getItem(2);
                end
            end
        
        % No parent means root has been reached and reachtube computation is complete
        elseif isempty(RTree.currNode.parent)
            % Update currNode now that tree search is complete
            RTree.currNode = treeNode.empty();
            
        % If a leaf node, check for left then right siblings, update currNode
        elseif RTree.currNode.parent.children.size > 1
            % If a left sibling exists
            lSib = RTree.currNode.parent.children.getItem(1);
            if lSib~=RTree.currNode
                if lSib.data.needsPost()==1
                    error('Left sibling should have Post computed already. Derp.');
                end
                % Combine reachsets computed in these sibling nodes
                RTree.currNode.data.reachUnion(lSib.data);
                % Pop left sibling
                RTree.currNode.parent.removeChild(1); % if lSib is referencing this node, then we might need to clear lSib or Matlab might not delete the object itself
                RTree.size = RTree.size-1;
            else
                RTree.currNode = RTree.currNode.parent.children.getItem(2); % A leaf node at this point should always be the left-most sibling; hopefully this is a handle to the correct node and not a copy; Try "RTree.currNode.parent.children.items{2}"
            end
        
        % If a leaf node and no siblings, update currNode to parent
        else
            % Combine reachsets with parent
            RTree.currNode.parent.data.reachUnion(RTree.currNode.data);
            % Update currNode to parent
            RTree.currNode = RTree.currNode.parent;
            % Delete child
            RTree.currNode.removeChild(1);
            RTree.size = RTree.size-1;
            RTree.height = RTree.height-1;
        end
    end
    % Update output
    cover.data.setSafeFlag(safetyEnum.RobustSafe);
    cover.data.reachtube = RTree.rootNode.data;
end