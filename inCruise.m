% Author:   Nicole Chan
% Created:  1/27/18
% Description: http://ctms.engin.umich.edu/CTMS/index.php?example=CruiseControl&section=SystemModeling
% SISO (1-dim) model with reference point x = vel = 10m/s
%
function funHand=inCruise()
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
%     m = 1000;   % [kg]
%     b = 50;     % [N.s/m]
%     A = -b/m;
%     B = 1/m;
    A =0.99501;
    B =0.0078125;

    % Closed-loop system
    ui = Fi*X+Gi;
    dX = A*X + B*ui;

end

function Jac=modelJacobian(X,Fi,Gi)
%     m = 1000;   % [kg]
%     b = 50;     % [N.s/m]
%     A = -b/m;
%     B = 1/m;
    A =0.99501;
    B =0.0078125;
    Jac = A+B*Fi;
end

function filename = MPCprob()
    filename = [];
end

function sol=MPCsol() 
    dim_x = 1;
    dim_u = 1;
    
    % Discrete time linear model in state-space rep: x_dot = Ax + Bu, y = Cx
    Ac =0.99501;
    Bc =0.0078125;
    dsys = c2d(ss(Ac,Bc,1,0),1); % Sampling period of 1s
    [A,B,~,~] = ssdata(dsys);
    
    % Cost function (infinity norm)
    Q = 1;
    R = 2;
    
    % Constraints
    XBND = 15;
    UBND = 1000;
    
    % Prediction horizon
    N = 10;
    
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
    expmpc = expmpc.simplify();
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
    xref = 8;
    xlow = 0;
    xup = 10;
    % Velocity should not be negative
    Uset.region(1,1) = Polyhedron('A',1,'b',xlow-xref);
    % Velocity should  not exceed 10m/s
    Uset.region(2,1) = Polyhedron('A',-1,'b',-(xup-xref));
    Uset.safe = zeros(2,1);
end

function X0=initialSet()
    xref = 8;
    x0 = 1;
    rad = 0.1;
    X0 = getBall(1,rad,x0-xref);
end

function vPar=verificationParams() 
    vPar.Thorizon = [0,10];      % verification time horizon
    vPar.simStep = 0.01;         % (fixed) simulation step size
    vPar.deltaStep = 1;          % sampling step size (same as what's used for discrete-time model and MPC period)
    vPar.maxPart = 6;            % maximum number of partitions taken per mode before returning an UNKNOWN result
    vPar.LipConst = 1;           % Lipschitz constant
    vPar.epsilonConst = 0;       % upper bound on simulation error
end