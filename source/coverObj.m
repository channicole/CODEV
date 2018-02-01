% Author:   Nicole Chan
% Created:  1/9/18
% Description: The 'coverObj' object is a subclass of 'reachtubeObj' that
% includes an extra parameter to indicate whether the entire reachtube has
% been computed and if so, if it is safe. The initial set of states for a
% 'coverObj' should always correspond to a partition of the user-defined
% inital set of states for verification purposes. The reachtube for a
% 'coverObj' should also cover the entire verification time horizon
% specified by the user. 
% In this way, 'coverObj' is distinguished from 'reachtubeObj', which is a 
% data structure that may be used for storing reachsets/tubes that only 
% cover some interval of the time horizon and use initial states that 
% result from subcomputations rather than the user-defined initial set
% Theta.
% 
% safeFlag: (1) reachtube computed and is safe; (0) reachtube computed and
% is not safe, so partitioning is needed; (-1) reachtube is not computed
% yet; (-2) exit verification process because the max number of partitions
% is reached without safe result
%
% NOTE: as of this iteration, we do not search for unsafe traces. 
%
classdef coverObj < handle
    properties
        reachtube@reachtubeObj  % reachtube data associated with this cover
    end
    properties (Access = protected)
%         reachtube@reachtubeObj  % reachtube data associated with this cover
        safeFlag@safetyEnum     % indicates whether a reachtube has been computed yet and if it is safe
    end
    
    methods
        % Constructor
        function cover = coverObj(Theta,x0,xi,MPCi,Reach,T,rad)
            if nargin ~= 7
                Theta = [];
                x0 = [];
                xi = [];
                MPCi = [];
                Reach = [];
                T = [];
                rad = [];
            end
            cover.reachtube = reachtubeObj(Theta,x0,xi,MPCi,Reach,T,rad);
            cover.safeFlag = safetyEnum.NeedReach;
        end
        
        % Set safeFlag
        function cover = setSafeFlag(cover,flag)
            if nargin ~= 2
                error('Incorrect number of arguments')
            end
            if cover.safeFlag ~= 1
                cover.safeFlag = flag;
            end
        end
        
        % Check safeFlag
        function flag = isSafe(cover)
            flag = cover.safeFlag;
        end
        
        % Return the superclass object
        function rtube = getReachtube(cover)
            rtube = cover.reachtube;
        end
   
        % Partition the current cover and returns the resulting 2^n covers
        % in a cell array
        function outCovers = partitionCover(inCover) 
            if nargin < 1
                error('Incorrect number of arguments')
            end
            
            [inx0,inrad,~] = inCover.reachtube.getProperties();
            n = length(inx0);
            outCovers = cell(2^n,1);
            
            if length(inx0)==n && length(inrad)==n 
                x0 = inx0;
                rad = inrad;
            else
                % Get ball approximation of the initial set
                [x0,rad] = poly2ball(inCover.reachtube.Theta);
            end
            
            if iscolumn(x0); x0 = x0'; end;
            if isrow(rad); rad = rad'; end;
            
            vals = [-0.5,0.5];
            X = cell(1, n);
            [X{:}] = ndgrid(vals);
            X = X(end : -1 : 1); 
            xcen = cat(n+1, X{:});
            xcen = reshape(xcen, [2^n, n]);
            xcen = xcen*diag(rad) + repmat(x0,2^n,1);
          
            rad = rad/2;
            for i=1:2^n
                % Get hyperbox that approximates the i-th partition of inCover.Theta
                boxTheta = Polyhedron('lb',xcen(i,:)'-rad,'ub',xcen(i,:)'+rad);
                % Get precise partition of inCover.reachtube.Theta using intersection
                Theta = inCover.reachtube.Theta & boxTheta;
                % Instantiate the coverObj for each partition
                outCovers{i} = coverObj(Theta,xcen(i,:)',[],inCover.reachtube.MPCi,boxTheta,inCover.reachtube.T,rad);
            end
        end
        
        % Concatenates computed cover from 'newCov' into 'cover'
        function cover = coverUnion(cover,newCov)
            if nargin~=2
                error('Incorrect number of arguments.')
            end
            cover.reachtube.reachUnion(newCov.reachtube);
        end
    end % end methods
end % end classdef