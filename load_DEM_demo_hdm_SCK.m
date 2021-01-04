function SCK=load_DEM_demo_hdm_SCK()
% This Demo is from SPM12
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_hdm_SCK.m 4804 2012-07-26 13:14:18Z karl $
  
% create model, generate data, set some parameters
% ==========================================================================
 
% set options
%--------------------------------------------------------------------------
clear M
M(1).E.linear = 0;                          % linear model
M(1).E.s      = 1;                          % smoothness
 
% level 1
%--------------------------------------------------------------------------
ip      = [1 2 5]';
pE      = spm_hdm_priors(1,5);              % parameters
np      = length(pE);
pC      = sparse(ip,ip,exp(-3),np,np);
pE(6)   = 0.02;
pE(7)   = 0.5;
 
M(1).n  = 4;
M(1).f  = 'spm_fx_hdm';
M(1).g  = 'spm_gx_hdm';
M(1).pE = pE;                               % prior expectation
M(1).pC = pC;                               % prior expectation
M(1).xP = exp(4);                           % prior expectation
M(1).V  = exp(8);                           % error precision
M(1).W  = exp(12);                          % error precision
M(1).ip = ip;
 
% level 2
%--------------------------------------------------------------------------
M(2).l  = 1;                                % inputs
M(2).V  = exp(0);                           % with shrinkage priors
 
% true parameters
%--------------------------------------------------------------------------
P       = M(1).pE;
P(ip)   = P(ip) - P(ip)/8;
 
 
% generate data
%==========================================================================
N         = 64;                             % length of data sequence
U         = exp(-([1:11] - 6).^2/(2.^2))/8; % this is the Gaussian cause
input     = zeros(1,N);
input([5 10 20 34 43 50]) = [1 0.8 1 0.2 .9 0.4];

U         = conv(U,input);
U         = U(1:N);

% U=Xin(N);
% tmp       = spm_DEM_generate(M,U,{P},{8,16},{16});
% U=full(tmp.Y);

DEM       = spm_DEM_generate(M,U,{P},{8,16},{16});

tmp=full(DEM.Y);
figure(10);
plot(tmp);

% t=0.1:0.1:6.4;
% for n=1:5
%  confound(n,:)=sin(n/4*pi*t);
% end
% confound=confound;
% rn=[1 -0.5 0.3 1.2 -0.1];
% 
% DEM.Y=DEM.Y+rn*confound;

spm_figure('GetWin','Simulated');
spm_DEM_qU(DEM.pU)
 
% Initialise SCK
%==========================================================================
DEM.M(1).E.nN = 32;
DEM.M(1).E.nN = 32;
SCK           = DEM;
 
% option to specify constrains on parameters (example, but not used here)
%-------------------------------------------------------------------------
cb(1,:) = [0.5,0.8];              % low and high bound
cb(2,:) = [0.3,0.5];              % low and high bound
cb(3,:) = [0.7,1.3];              % low and high bound
cb(4,:) = [0.29,0.35];            % low and high bound
cb(5,:) = [0.32,.4];              % low and high bound
cb(6,:) = [0.01,0.04];            % low and high bound
cb(7,:) = [0.5,0.6];              % low and high bound
 
SCK.M(1).ip = ip;  % indices of model parameters to be estimated
SCK.M(1).cb = [];  % option to specify constrain on parameters values [min max]
SCK.M(2).v  = 0;   % input initial condition
SCK.M(2).V  = 5;   % input noise precision (fixed)
SCK.M(1).xP = eye(4)*1e-1^2;      % state error covariance matrix
SCK.M(1).uP = eye(1)*1e-2^2;      % input error covariance matrix
SCK.M(1).wP = eye(np)*1e-10^2;    % parameter error covariance matrix
SCK.M(1).pC = eye(np)*1e-10^4;    % parameter error covariance matrix
SCK.M(1).f  = 'spm_fx_hdm_sck';   % state equations rewritten for matrix operations
SCK.M(1).g  = 'spm_gx_hdm_sck';   % observation equations rewritten for matrix operations
SCK.M(1).Q  = {}; 
SCK.M(1).E.nN    = 32;            % max number of iteration of SCKF-SCK algorithm
SCK.M(1).E.nD    = 4;             % time-steps
SCK.M(1).E.Itol  = 1e-3;          % convergence tolerance value for SCKF_SCK algorithm
SCK.M(1).E.RM    = [5e3 1e6];     % scaling parameter for Robbins-Monro approximation 

% save SCK SCK;