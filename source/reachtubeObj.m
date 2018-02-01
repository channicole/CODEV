% Author:   Nicole Chan
% Created:  10/6/17
% Description: The 'reachtubeObj' object represents a polyhedral subset of 
% the state-space that contains all the reachable states over the time 
% horizon 'T'.
% It contains all the reachsets R_i, for all i in T, where 
% R_i is the reachable states over t in [t_{i-1},t_{i}).
%
classdef reachtubeObj < handle
    properties (Access = public)
        Theta   % Polyhedron: initial set of states
%         x0      % n-length vector: initial state in Theta to simulate from
        xi      % k-length array of n-length vectors: nominal trajectory simulated
        MPCi    % k-length array of MPC-objects: affine control solution applied
%         Reach   % k-length array of Polyhedron: reachsets (Reach(T(1))=Theta)
        T       % k-length vector: discrete time steps
%         rad     % n-length vector: radius of the ball approximation of Theta
    end
    % The following members represent the hyperrectangular approximation of
    % Theta. The Polyhedron representation is stored in Reach(1,:). Its
    % center state and radius is stored in x0 and rad, respectively.
    properties (Access = protected)
        x0      % n-length vector: initial state in Theta to simulate from
        Reach   % k-length array of Polyhedron: reachsets (Reach(T(1))=Theta)
        rad     % n-length vector: radius of the ball approximation of Theta
    end
    
    methods
        function reach = reachtubeObj(Theta,x0,xi,MPCi,Reach,T,rad) % constructor
            if nargin == 7
                if isrow(x0); x0 = x0'; end;
                if iscolumn(xi); xi = xi'; end;
                if isrow(rad); rad = rad'; end;
                reach.Theta = Theta;
                reach.x0 = x0;
                reach.xi = xi;
                reach.MPCi = MPCi;
                reach.Reach = Reach;
                reach.T = T;
                reach.rad = rad;
            else
                error('Incorrect number of arguments')
            end
        end
        % Concatenates reachtube from 'newReach' into 'reach'
        function reach = reachUnion(reach,newReach)
            if nargin~=2
                error('Incorrect number of arguments.')
            end
            % Copy time horizons and truncate to appropriate length
            % corresponding with reachsets
            T1 = reach.T(1:length(reach.Reach));
            T2 = newReach.T(1:length(newReach.Reach));
            [~,i,j] = intersect(T1,T2);
%             [~,i,j] = intersect(reach.T,newReach.T);
            for k=1:length(i)
                % Convex hull method:
                rUnion = PolyUnion([reach.Reach(i(k),:);newReach.Reach(j(k),:)]);
                try
                    reach.Reach(i(k),:) = rUnion.convexHull;
                catch
                    error('ConvexHull operation throwing errors.');
                end
            end
            k = setdiff((1:length(newReach.T))',j);
            if isempty(reach.Reach)
                reach.Reach = newReach.Reach(k,:); % this should only happen for coverUnion when partitioning for different modes initially
            else
                reach.Reach(end+1:end+length(k),:) = newReach.Reach(k,:);
            end
        end
        % Updates the reachtube/set members
        function reach = updateReach(reach,newReach)
            if size(newReach,1) == 1 
                [reach.x0,reach.rad,reach.Reach,~] = poly2ball(newReach);
                reach.rad = reach.rad';
            elseif ~isempty(newReach) % if newReach has reachsets for multiple time-steps
                reach.Reach = newReach;
                [reach.x0,reach.rad,~,~] = poly2ball(newReach(1));
                reach.rad = reach.rad';
            else
                warning('Input is a handle to an empty object, so member not updated.')
            end
        end
        % Returns a bool on whether the reachtube/set has been computed yet
        function flag = needsPost(inReach)
            if size(inReach.Reach,1)<2
                flag = 1; % based on the assumption that min length of inReach.T is 2, so if computed already then there should be at least two corresponding sets
            elseif size(inReach.Reach,1)<length(inReach.T)
                flag = 2;
            else
                flag = 0;
            end
        end
        % Returns the protected properties
        function [x0,rad,Reach] = getProperties(inReach)
            x0 = inReach.x0;
            rad = inReach.rad;
            Reach = inReach.Reach;
        end
    end % end methods
end % end classdef