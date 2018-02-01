% Name:     MPC_to_mpLP.m
% Author:   Nicole Chan
% Date created: 10/2/17
% Description: Recasts MPC problem into mp-LP and stores it in a Opt object
% from MPT toolbox. Additionally solves the problem and returns solution.
% This part was originally defined in MPCsol_mpLP_solve.m.
% It should be noted that for a m-dimensional input, the resulting input
% vector is represented by a m*Nu-dimensional vector:
% [u1(0),u2(0),...,um(0),u1(1),u2(1),..,um(1),u1(Nu-1),...,um(Nu-1)]
% 
function [mplp_prob,mplp_sol,Fi,Gi]=MPC_to_mpLP(modelParams)
%% Setup
% Load variables from the file whose name is the string modelParams
try
    load(modelParams,'dim_x','dim_u','dim_y','Nu','Ny','P','Q','R',...
        'xbnd','ubnd','ybnd','A','B','C');
catch
%     error('Could not load parameters for MPC problem. Please follow template for MPCprob() function.');
    % Struct containing the parameters
    try
        A =     modelParams.A;
        B =     modelParams.B;
        C =     modelParams.C;
        dim_x = modelParams.dim_x;
        dim_u = modelParams.dim_u;
        dim_y = modelParams.dim_y;
%         Nc =    modelParams.Nc;
        Nu =    modelParams.Nu;
        Ny =    modelParams.Ny;
        P =     modelParams.P;
        Q =     modelParams.Q;
        R =     modelParams.R;
        xbnd =  modelParams.xbnd;
        ubnd =  modelParams.ubnd;
        ybnd =  modelParams.ybnd;
    catch
        error('Please specify input file or the appropriate parameters in a struct.');
    end
end


dim_e = Ny+Nu;              % number of slack variables
dim_z = dim_e + dim_u*Nu;   % dimension of the recasted mpLP state vector

% Process constraints
if size(xbnd,1)==1 && size(xbnd,2)==2
    xmin = xbnd(1)*ones(dim_x,1);
    xmax = xbnd(2)*ones(dim_x,1);
elseif size(xbnd,1)==dim_x && size(xbnd,2)==2
    xmin = xbnd(:,1);
    xmax = xbnd(:,2);
elseif ~isempty(xbnd)
    error('Incorrect format for state constraints');
end

if size(ubnd,1)==1 && size(ubnd,2)==2
    umin = ubnd(1)*ones(dim_u,1);
    umax = ubnd(2)*ones(dim_u,1);
elseif size(ubnd,1)==dim_u && size(ubnd,2)==2
    umin = ubnd(:,1);
    umax = ubnd(:,2);
elseif ~isempty(ubnd)
    error('Incorrect format for input constraints');
end

if size(ybnd,1)==1 && size(ybnd,2)==2
    ymin = ybnd(1)*ones(dim_y,1);
    ymax = ybnd(2)*ones(dim_y,1);
elseif size(ybnd,1)==dim_y && size(ybnd,2)==2
    ymin = ybnd(:,1);
    ymax = ybnd(:,2);
elseif ~isempty(ybnd)
    error('Incorrect format for output constraints');
end

%% Optimization program def: min J(z) = c*z, such that G*z <= w + F*x
% Objective function for recasted problem
c = ones(1,dim_z);
c((dim_e+1):end) = zeros(1,dim_u*Nu);

% % % Constraints: % % %
% Slack_x
numRows = 2*dim_x; 
G1 = zeros(Ny*numRows,dim_z);
w1 = zeros(Ny*numRows,1);
F1 = zeros(Ny*numRows,dim_x);
for i=1:Ny
    G1(1+(i-1)*numRows:i*numRows,i) = -ones(numRows,1); % slack variables
    F1(1+(i-1)*numRows:i*numRows,:) = [Q*A^i; -Q*A^i];  % F1*x(t)
    % -Q*A^i*B*[u1,...,uN]'
    for j=i:i+numRows/2
        for k=0:i-1
            G1(1+(i-1)*numRows:i*numRows,dim_e+i-k:dim_e+i-k+(dim_u-1)) = [-Q*A^k*B; Q*A^k*B];
        end
    end
end

% Slack_u
numRows = 2*dim_u;
G2 = zeros(Nu*numRows,dim_z);
w2 = zeros(Nu*numRows,1);
F2 = zeros(Nu*numRows,dim_x);
for i=1:Nu
    G2(1+(i-1)*numRows:i*numRows,i) = -ones(numRows,1); % slack variables
    G2(1+(i-1)*numRows:i*numRows,dim_e+i:dim_e+i+(dim_u-1)) = [-R; R]; % -R*[u1,...,uN]
end

% State contraints
if isempty(xbnd)
    G3 = [];
    w3 = [];
    F3 = [];
else
    numRows = 2*dim_x;
    G3 = zeros(Ny*numRows,dim_z);
    w3 = zeros(Ny*numRows,1);
    F3 = zeros(Ny*numRows,dim_x);
    for i=1:Ny
        w3(1+(i-1)*numRows:i*numRows) = [-xmin;xmax];
        F3(1+(i-1)*numRows:i*numRows,:) = [A^i;-A^i];
        for j=0:i-1
            G3(1+(i-1)*numRows:i*numRows,dim_e+(i-j-1)*dim_u+1:dim_e+(i-j)*dim_u) = [-A^j*B;A^j*B];
        end
    end
end

% Output constraints 
if isempty(ybnd)
    G4 = [];
    w4 = [];
    F4 = [];
else
    numRows = 2*dim_y;
    G4 = zeros(Ny*numRows,dim_z);
    w4 = zeros(Ny*numRows,1);
    F4 = zeros(Ny*numRows,dim_x);
    for i=1:Ny
        w4(1+(i-1)*numRows:i*numRows) = [-ymin;ymax];
        F4(1+(i-1)*numRows:i*numRows,:) = [C*A^i;-C*A^i];
        for j=0:i-1
            G4(1+(i-1)*numRows:i*numRows,dim_e+i-j:dim_e+i-j+(dim_u-1)) = [-C*A^j*B;C*A^j*B];
        end
    end
end

% Input constraints
if isempty(ubnd)
    G5 = [];
    w5 = [];
    F5 = [];
else
    numRows = 2*dim_u;
    G5 = zeros(Nu*numRows,dim_z);
    w5 = zeros(Nu*numRows,1);
    F5 = zeros(Nu*numRows,dim_x);
    for i=1:Nu
        w5(1+(i-1)*numRows:i*numRows) = [-umin;umax];
        G5(1+(i-1)*numRows:i*numRows,dim_e+(i-1)*dim_u+1:dim_e+i*dim_u) = [-eye(dim_u);eye(dim_u)];
    end
end

% Collect into general form (with nonzero slack constraints)
G = [G1;G2;G3;G4;G5;-eye(dim_e),zeros(dim_e,dim_u*Nu)];
w = [w1;w2;w3;w4;w5;zeros(dim_e,1)];
F = [F1;F2;F3;F4;F5;zeros(dim_e,dim_x)];


%% Solve mpLP using mpt toolbox
% Formulate optimization problem using Opt class using convention: 
% min f'*x s.t. A*x<=b+pB*theta
mplp_prob = Opt('f',c','A',G,'b',w,'pB',F);
mplp_sol = mplp_prob.solve();

% If only using mpLP solution for MPC, remove slack vars from solution
Fi = cell(1,mplp_sol.xopt.Num); 
Gi = cell(1,mplp_sol.xopt.Num);
if mplp_sol.xopt.Num>0
    optimizer = mplp_sol.xopt.Set.getFunction('primal');
    for i = 1:length(optimizer)
        Fi{i} = optimizer(i).F(Ny+Nu+1:Ny+Nu+dim_u,:);
        Gi{i} = optimizer(i).g(Ny+Nu+1:Ny+Nu+dim_u);
    end
end
