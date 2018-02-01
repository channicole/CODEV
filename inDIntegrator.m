% Author:   Nicole Chan
% Created:  1/23/18
% Description: This is the template user-defined input file
%
function funHand=inDIntegrator()
    % Create function handles to the hybrid model components, MPC
    % formulation, and verification parameters
    funHand.flowEq = @flowModel;        % ODEs in each discrete mode (possibly nonlinear)
    funHand.modelJac = @modelJacobian;  % Jacobian in each mode
    funHand.unsafeStates = @unsafeSetWrapper;   % Sets of unsafe states
    funHand.initStates = initialSet;    % Initial set of states for verification
    funHand.vPar = verificationParams;  % Time horizon, error terms, sim step-size
    funHand.MPCprob = MPCprob;          % Linear, discrete-time plant model used to design MPC sol
    funHand.MPCsol = MPCsol;
end

function dX=flowModel(t,X,Fi,Gi)
    % Continuous time model corresponding with the discrete-time given in MPCprob
    A = [0,1;0,0];
    B = [0;1];

    % Closed-loop system
    ui = Fi*X+Gi;
    dX = A*X + B*ui;

end

function Jac=modelJacobian(X,Fi,Gi)
    A = [0,1;0,0];
    B = [0;1];
    Jac = A+B*Fi;
end

function filename = MPCprob()
    filename = [];
end

function sol=MPCsol()
    dim_x = 2;
    dim_u = 1;
    
    % Discrete time linear model in state-space rep: x_dot = Ax + Bu, y = Cx
    A = [1,1;0,1];
    B = [0.5;1];
    
    % Cost function (infinity norm)
    Q = [1,1;0,1];
    R = 0.8;
    
    % Constraints
    XBND = 10;
    UBND = 1;
    
    % Prediction horizon
    N = 2;
    
    % Use MPT libraries to obtain explicit MPC solution and then copy into
    % our custom struct
    model = LTISystem('A',A,'B',B);
    model.x.min = -XBND*ones(dim_x,1);
    model.x.max = XBND*ones(dim_x,1);
    model.u.min = -UBND*ones(dim_u,1);
    model.u.max = UBND*ones(dim_u,1);
    model.x.penalty = InfNormFunction(Q);
    model.u.penalty = InfNormFunction(R);
    mpc = MPCController(model,N);
    expmpc = mpc.toExplicit();
    optimizer = expmpc.optimizer;
    
    % Custom struct:
    sol.Pn = optimizer.Set;
    sol.Fi = {}; 
    sol.Gi = {};
    for i = 1:optimizer.Num
        sol.Fi{i} = optimizer.Set(i).Functions('primal').F(1:dim_u,:);
        sol.Gi{i} = optimizer.Set(i).Functions('primal').g(1:dim_u,:);
    end
end

% Use a wrapper function to return both the cell array of actual unsafe
% sets and a function handle to instantiate these sets when there is a
% constraint on non-explicit control input variable OR whenever the unsafe
% region is parameterized based on the discrete mode invariants
function unsafeArray = unsafeSetWrapper(numRegions)
    unsafeArray.array{numRegions,1} = [];
    unsafeArray.funHand = @unsafeSet;
end
function Uset=unsafeSet(F,G)
    % Currently unsafe sets are represented by Polyhedron objects
    rad = [10,10];
    Uset.region(1,1) = getBall(2,rad); % state constraints
    Uset.region(2,1) = Polyhedron('A',[F;-F],'b',[1-G;1+G]); % input constraints
    Uset.safe = ones(2,1);
end

function X0=initialSet()
    x0 = [9,-5];
    rad = [0.01,0.01];
    X0 = getBall(2,rad,x0);
end

function vPar=verificationParams() 
    vPar.Thorizon = [0,9];       % verification time horizon
    vPar.simStep = 0.01;         % (fixed) simulation step size
    vPar.deltaStep = 1;          % sampling step size (same as what's used for discrete-time model and MPC period)
    vPar.maxPart = 6;            % maximum number of partitions taken per mode before returning an UNKNOWN result
    vPar.LipConst = 1;           % Lipschitz constant
    vPar.epsilonConst = 0;       % upper bound on simulation error
    
end