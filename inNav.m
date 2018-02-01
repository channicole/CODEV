% Author:   Nicole Chan
% Created:  1/27/18
% Description: 
%
function funHand=inNav()
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
    A = [1,0,0.1,0;0,1,0,0.1;0,0,0.8870,0.0089;0,0,0.0089,0.8870];
    B = [0,0;0,0;1,0;0,1];

    % Closed-loop system
    ui = Fi*X+Gi;
    dX = A*X + B*ui;
end

function Jac=modelJacobian(X,Fi,Gi)
    A = [1,0,0.1,0;0,1,0,0.1;0,0,0.8870,0.0089;0,0,0.0089,0.8870];
    B = [0,0;0,0;1,0;0,1];
    Jac = A+B*Fi;
end

function filename = MPCprob()
    filename = [];
end

function sol=MPCsol() 
    dim_x = 4;
    dim_u = 2;
    
    % Discrete time linear model in state-space rep: x_dot = Ax + Bu, y = Cx
    Ac = [1,0,0.1,0;0,1,0,0.1;0,0,0.8870,0.0089;0,0,0.0089,0.8870];
    Bc = [0,0;0,0;1,0;0,1];
    dsys = c2d(ss(Ac,Bc,eye(4),0),1); % Sampling period of 1s
    [A,B,~,~] = ssdata(dsys);
    
    % Cost function (infinity norm)
    Q = eye(4);
    R = eye(2);
    
    % Constraints
    XBND = 40;
    UBND = 100;
    
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
%     load('Nav_expmpc.mat');
    mpc = MPCController(model,N);
    expmpc = mpc.toExplicit();
%     expmpc = expmpc.simplify();
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

function Uset=unsafeSet(~,~)
    xref = [10;12;0;0];
    Uset.region(1,1) = getBall(4,[0.15,0.8,0,0],[5.65;5.8;0;0]-xref);
    Uset.region(2,1) = getBall(4,[1.65,2.3,0,0],[5.15;6.3;0;0]-xref);
    Uset.region(3,1) = getBall(4,[1.5,1.5,0,0],[3.5;1.5;0;0]-xref);
    Uset.safe = zeros(3,1);
end

function X0=initialSet()
    xref = [10;12;0;0];
    x0 = [8;10;1;1];
    rad = 0.1;
    X0 = getBall(4,rad,x0-xref);
end

function vPar=verificationParams() 
    vPar.Thorizon = [0,2];       % verification time horizon
    vPar.simStep = 0.01;         % (fixed) simulation step size
    vPar.deltaStep = 1;          % sampling step size (same as what's used for discrete-time model and MPC period)
    vPar.maxPart = 6;            % maximum number of partitions taken per mode before returning an UNKNOWN result
    vPar.LipConst = 1;           % Lipschitz constant
    vPar.epsilonConst = 0;       % upper bound on simulation error
    
end