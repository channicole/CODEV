% Author:   Nicole Chan
% Created:  1/27/18
% Description: 
%
function funHand=inMagPointerC()
    % Create function handles to the hybrid model components, MPC
    % formulation, and verification parameters
    funHand.flowEq = @flowModel;        % ODEs in each discrete mode (possibly nonlinear)
    funHand.modelJac = @modelJacobian;  % Jacobian in each mode
    funHand.unsafeStates = unsafeSet;   % Sets of unsafe states
    funHand.initStates = initialSet;    % Initial set of states for verification
    funHand.vPar = verificationParams;  % Time horizon, error terms, sim step-size
    funHand.MPCprob = MPCprob;          % Linear, discrete-time plant model used to design MPC sol
    funHand.MPCsol = MPCsol;
end

function dX=flowModel(t,X,Fi,Gi)
    A = [2.6629,-1.1644,0.66598;2,0,0;0,0.5,0];
	B = [0.25;0;0];
    
    % Closed-loop system
    ui = Fi*X+Gi;   % scalar input
    dX = A*X + B*ui;
end

function Jac=modelJacobian(X,Fi,Gi)
    A = [2.6629,-1.1644,0.66598;2,0,0;0,0.5,0];
	B = [0.25;0;0];
    Jac = A+B*Fi;
end

function filename = MPCprob()
    filename = [];
end

function sol=MPCsol() 
    dim_x = 3;
    dim_u = 1;
    
    % Linearized continuous-time model
    Ac = [2.6629,-1.1644,0.66598;2,0,0;0,0.5,0];
	Bc = [0.25;0;0];

    % Discrete time linear model in state-space rep: x(k+1) = Ax(k) + Bu(k)
    dsys = c2d(ss(Ac,Bc,eye(dim_x),0),1); 
    [A,B,~,~] = ssdata(dsys);
    
    % Cost function (infinity norm)
    Q = eye(dim_x);
    R = eye(dim_u);
    
    % Constraints
    XBND = 5;
    UBND = 40;
    
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

function Uset=unsafeSet(F,G)
    rad = 5;
    Uset.region(1) = getBall(3,rad,zeros(3,1));
    Uset.safe = 1;
end

function X0=initialSet()
    xref = [1;2;2];
    x0 = zeros(3,1);
    rad = 0.1*ones(3,1);
    X0 = getBall(3,rad,x0-xref);
end

function vPar=verificationParams() 
    vPar.Thorizon = [0,1.1];     % verification time horizon
    vPar.simStep = 0.01;         % (fixed) simulation step size
    vPar.deltaStep = 1;          % sampling step size (same as what's used for discrete-time model and MPC period)
    vPar.maxPart = 6;            % maximum number of partitions taken per mode before returning an UNKNOWN result
    vPar.LipConst = 1;           % Lipschitz constant
    vPar.epsilonConst = 0;       % upper bound on simulation error
    
end