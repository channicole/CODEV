% Author:   Nicole Chan
% Created:  10/9/17
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a continuous (single-mode) system:
% 1. Computes  a 'center' state to simulate from
% 2. Simulate for length of time T (T should be a vector of all simulation
% time steps)
% 3. Compute discrepancy function and bloat the simulation trace
%
function [Reachtube,safeflag]=computePost(Theta,T,modelFun,modelMPC,MPCi,modelJac,unsafeSet,epsilonConst,LipConst,checkSafe)
    %% Setup
    if nargin == 6
        epsilonConst = 0;
        LipConst = 0;
    elseif nargin == 7
        LipConst = 0;
    end

    % Initialize safeflag to unknown
    safeflag = -1;
    
    % Get ball approximation of the initial set
    [x0,deltaConst,ball0] = poly2ball(Theta,modelMPC.Pn,MPCi.ind);
    Reachtube = reachtubeObj(Theta,x0,x0',MPCi,ball0,T,deltaConst);

    %% Simulate
    [sim.T,sim.X] = ode45(@(t,x) modelFun(t,x,Reachtube.MPCi.Fi,Reachtube.MPCi.Gi),T,x0); 

    %% Compute discrepancy (ComputeLDF algorithm from ATVA15)
    x_dim = length(x0);
    b = zeros(length(sim.T),1);
    b(1) = max(deltaConst);
    beta = b;
    Delta = b(1);
    dt = sim.T(2)-sim.T(1);
    % % TODO: add parameter to input file for LTI flag
    % If the Jacobian is static (e.g. LTI model), eliminate unnecessary function calls
    if ~isa(modelJac,'function_handle')
        evalJac = 0;
        J = modelJac;
        lambda = max(eig(J+J')/2);
    else
        evalJac = 1;
    end
    for i=2:length(sim.T)
        if evalJac
            dia = (Delta + epsilonConst)*exp(LipConst*dt);% amount of coarse bloating
            S = Polyhedron([sim.X(i-1,:);sim.X(i,:)]);  % coarse overapprox of Reach over dt
            if dia ~= 0
                S = plus(S,getBall(x_dim,dia));
            end
            [centerS,~,~] = poly2ball(S);               % get center state of S
            J = modelJac(centerS,Reachtube.MPCi.Fi,Reachtube.MPCi.Gi); % get Jacobian
            lambda = max(eig(J+J')/2);                  % 
        end
        err = 0;                                    % ignore line 8 for now
        b(i) = lambda + err/2;                      % store local discrepancy
        Delta = (Delta + epsilonConst)*exp(b(i)*dt);% update error
        beta(i) = beta(i-1)*exp(b(i)*dt);           % missing from Alg 2
    end
    
    %% Bloat simulation trace
    rtube(length(sim.T),1) = Polyhedron();
    for i=1:length(sim.T)
        rtube(i) = getBall(x_dim,beta(i),sim.X(i,:));
    end

    %% Check safety
    if checkSafe
        [safeflag,ind] = checkSafety(rtube,unsafeSet);
        if ~safeflag
            disp(['Unsafe at t=',num2str(T(ind))]);
            return
        end
    end
    
    %% Update Reachtube parameters with computed reachtube data
    Reachtube.updateReach(rtube);
    Reachtube.MPCi = repmat(Reachtube.MPCi,length(sim.T),1);
    Reachtube.xi = sim.X;
end