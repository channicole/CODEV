% Name:     MPCsolve.m
% Author:   Nicole Chan
% Date created: 10/3/17
% Description: Loads car models for different parameters, solves MPC
% solution, 
%
% UPDATE 1/14/18: udpated Fi, Gi to only output relevant rows

function [Pn,Fi,Gi,activeConstraints,Phard,details]=MPCsolve(MPCprob)
% Solve for MPC solution
[~,mplp_sol,Fi,Gi]=MPC_to_mpLP(MPCprob);

Pn = mplp_sol.xopt.Set;
% ACTIVECONSTRAINTS ONLY COMPUTED USING 'MPLP' SOLVER NOT 'PLCP'
try 
    activeConstraints = mplp_sol.mplpsol.activeConstraints;
catch
    activeConstraints = {};
end
Phard = mplp_sol.xopt.Domain;


% Cost
details.Ai = cell(1,mplp_sol.xopt.Num);
details.Bi = cell(1,mplp_sol.xopt.Num);
details.Ci = cell(1,mplp_sol.xopt.Num);
if mplp_sol.xopt.Num>0
    cost = mplp_sol.xopt.Set.getFunction('obj');
    for i = 1:length(cost)
        if isa(cost(i), 'QuadFunction')
            details.Ai{i} = cost(i).H;
        end
        details.Bi{i} = cost(i).F;
        details.Ci{i} = cost(i).g;
    end
end

end