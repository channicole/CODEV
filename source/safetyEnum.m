% Author:   Nicole Chan
% Created:  1/9/18
% Description: 
%
classdef safetyEnum < int8
    enumeration
        RobustSafe  (1)
        NotSafe     (0)
        NeedReach   (-1)
        ExitPartitionBnd (-2)
        ExitInfeasible (-3)
    end
end