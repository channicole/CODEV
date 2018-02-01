% Author:   Nicole Chan
% Created:  11/10/17
% Description: Checks the intersection between the reachtube (a 
% reachtubeObj object) and unsafe sets (stored in a function that unsafeSet
% points to). Returns true if the intersection is empty, and false 
% otherwise.
%
% NOTE: both reachtube and unsafeSet should be arrays of Polyhedron()
%
function [safeflag,ind]=checkSafety(reachtube,unsafeSet)
    safeflag = 1;
    ind = -1;
    % Check each property
    for i = 1:length(unsafeSet.safe)
        % If the property specified is the safe set (i.e. complement of
        % unsafe), then check the reachtube is contained in this set
        if unsafeSet.safe(i)
            for j = 1:length(reachtube)
                if ~isEmptySet(reachtube(j)\unsafeSet.region(i,:))
%                 if ~unsafeSet.region(i).contains(reachtube(j))
                    safeflag = 0;
                    ind = j;
                    return;
                end
            end
        % Else if the property specified is the unsafe set, check the 
        % reachtube does not intersect with this set
        else
            for j = 1:length(reachtube)
                if any(~isEmptySet(unsafeSet.region(i,:) & reachtube(j)))
                    safeflag = 0;
                    ind = j;
                    return;
                end
            end
        end  
    end
end