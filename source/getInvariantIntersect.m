% Author:   Nicole Chan
% Created:  1/20/18
% Description: Checks the intersection between the reachtube (a 
% reachtubeObj object) and unsafe sets (stored in a function that unsafeSet
% points to). Returns true if the intersection is empty, and false 
% otherwise.
%
% NOTE: both reachtube and unsafeSet should be arrays of Polyhedron()
%
% RETURNS: 'outRegions' is the reachtube (Polyhedron-vector) corresponding
% with 'currRegion'
% 'nextRegions' is an array of reachsets (Polyhedron) that correspond with
% a transition
% 'nextInd' is an array of indices in 'outRegions' that correspond with the
% set in 'nextRegions' (temporal)
% 'regInd' is an array of indices in 'invRegions' (spatial)
%
function [outRegions,nextRegions,nextInd,regInd]=getInvariantIntersect(reachtube,invRegions,currRegion)
    % Setup
    if nargin==3
        % Check the currRegion is valid
        test = invRegions.contains(currRegion);
        if length(find(test))~=1
            error('Check correctness of currRegion input arg.');
        end
    elseif nargin==2
        % Check currRegion can be uniquely obtained from the invariant space
        currRegion = invRegions & reachtube(1);
        if length(find(~isEmptySet(currRegion))) > 1
            error('This initial set covers more than one region.');
        end
        currRegion = currRegion(~isEmptySet(currRegion));
    else
        error('Incorrect input arguments.');
    end
    
    % Find intersection with current mode (used for reachtube computation)
    intReach = reachtube & currRegion;
    
    % ind = index of first empty intersection-->a transition MUST be taken
    ind = find(isEmptySet(intReach),1);
    outRegions = intReach(1:ind)';
    
    % Get intersections with other modes up to reachtube(ind)
    j = 1;
    for i=2:ind
        intRegion = invRegions & reachtube(i);
        k = find(~intRegion.isEmptySet & intRegion.isFullDim);
        if ~isempty(k)
            nextRegions(j,:) = intRegion(k);
            nextInd(j) = i;
            regInd(j) = k(1);
            j = j+1;
        end
    end
end