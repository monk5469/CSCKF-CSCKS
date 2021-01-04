% DEMO for deconvolution using SCKS algorithm
% Created by Martin Havlicek 2011, Albuquerque
% Improved by Haifeng Wu 2019, Yunnan Minzu University
% The improved Demo is for confounds and the details can be the follosing
% paper.

% [1] Haifeng Wu, Mingzhi Lu, Yu Zeng. State Estimation of Hemodynamic 
%     Model for fMRI Under Confounds: SSM Method. IEEE Journal of 
%     Biomedical and Health Informatics, Volume 24 , Issue 3 , March 2020.
%==================================================================

clear all; 
close all;

%%%%%%%%%%%%%% I. This part is Martin's Demo, where the parameters %%%%%%%
%%%%%%%%%%%%%% in DCM is set. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Here add path to spm12   !!!!!!!!
% addpath(genpath('....\spm12\'));

% 2. Load data:
% data (=time course) should be for example the first eigen-variate from the
% ROI. Data should be already corrected for slow fluctuations and motion.
% We further assume data having zero-mean (or even better adjusted to zero baseline in rest)
load data.mat
DCM.Y = data;
 
T   = length(DCM.Y);                             % number of time point

% 3. Here write the TR of your data:
TR  = 2;                                         % sampling period

% scale data to have particular peak to peak range e.g. 2 % signal
% change
scale = peak2peak(DCM.Y');

pk2pk = 2;   
if scale > pk2pk
    scale   = pk2pk/max(scale,pk2pk);
    DCM.Y   = DCM.Y*scale;
else
    scale   = pk2pk/min(scale,pk2pk);
    DCM.Y   = DCM.Y*scale;
end

% Specify generative model for inversion (DCM)
% =========================================================================
% set inversion parameters
% -------------------------------------------------------------------------

DCM.M(2).v    = 0;
DCM.M(1).E.TR = TR; 
DCM.M(1).E.dt = 1; 
n             = 1;
% Specify hyper-priors on precisions
% -------------------------------------------------------------------------
W             = exp(spm_vec(sparse(1:n,1,(6 - 16),n,5) + 16));
DCM.M(1).xP   = exp(6);
DCM.M(1).V    = 1;             % prior log precision (noise)
DCM.M(1).W    = diag(W);       % fixed precision (hidden-state)
DCM.M(2).V    = exp(8);        % fixed precision (hidden-cause) (not used)


% Full connectivity inversion
%==========================================================================
F  = ones(n,n);
B  = zeros(n,n,0);
C  = zeros(n,1);
D  = zeros(n,n,0);
options.endogenous = 0;  
[pE pC ill pW]     = spm_dcm_fmri_priorsM(F,B,C,D,options);   % generating parameter priors for model inversion
                                                              % 
                                                              
DCM.M(1).pE = pE;
DCM.M(1).pC = pC;
DCM.M(1).cb = [];
DCM.M(1).Nc = n^2;

DCM.M(1).W  = diag(W)/TR*1e3;        %!!! precision on neuronal (and hemodynamic) state noise - IMPORTANT
                                   % might need to be adjustment to your (it defines smoothness of the estimate)
                                  
DCM.M(2).V  = eye(n);
DCM.M(1).V  = eye(n)*DCM.M(1).V;
DCM.M(2).v  = 0;
DCM.M(1).x  = zeros(n,4);

% defining paramters (their covariances) for estimation:
nA   = ones(n);                                  % full connectivity matrix
nA   = diag(nA(:));
hP   = diag(repmat([0 1 0 1 0]',n,1));           % select hemodynamic paremeters for estimation;
DCM.M(1).pC = pC.*blkdiag(nA,eye(n),hP,1);       % prior covariance on parameters

DCM.M(1).pE = pE;
DCM.M(1).pP = pW*1e-6/3;  % 

DCM.M(1).wP = pC*1e-2;
DCM.M(1).uP = 0; 
DCM.M(1).xP = blkdiag(eye(4*n)*1e-1)/3; %

pE          = spm_vec(pE);
ip          = [1:length(pE)]; 

DCM.M(1).ip = ip;
DCM.M(1).l  = n;
DCM.M(1).Q  = {speye(DCM.M(1).l,DCM.M(1).l)};   % if Q is specified then algorithm performs
                                                % estimation of measurement noise covariance 
%DCM.M(1).Q  = [];                              % if presion on measurement noise
                                                % is known then Q = [];
DCM.M(1).E.nN    = 20;                   % max number of iteration of SCKF-SCKS algorithm
DCM.M(1).E.Itol  = 1e-5;                 % convergence tolerance value for SCKF_SCKS algorithm
DCM.M(1).E.nD    = 5;                   %numeber of integration steps.!!! set to match the TR but it can be also higher number (but it might need other adjustments)
DCM.M(1).E.RM    = [1e2 1e8];
DCM.M(1).VB.N    = 5;                    % number of VB iteration during one SCKF-SCKS run
DCM.M(1).VB.l    = 1 - exp(-6);          % scaling parameter for VB algorithm, 
                                         % controls dynamics

fx = ['@(x,v,P) spm_fx_fmri4(x,v,P)*',num2str(TR)]; % functions of hemodynamic model for model inversion by SCKS 
gx = ['@(x,v,P) spm_gx_fmri4(x,v,P)'];

DCM.M(1).f    = eval(fx);
DCM.M(1).g    = eval(gx);
DCM.pU        = [];
DCM.pU.x{1}   = [];
%==================================================================



%%%%%%%%%%%%%% II. This part is to simulate confounds %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%% generate confounds,old edition %%%%%%%%%%%%%%%%%
% tmp=0.1:0.1:6.4;
% for ntmp=1:6
%     confounds(ntmp,:)=sin(ntmp/4*pi*tmp);
% end
% rn=[1.2 -0.2 0.5 -0.4 0.3 0.1];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%%%%%%%%%  confounds, 20190314 edition %%%%%%%%%%%%%%%%%
N=64;
rn=[1.2 -0.2 0.5 -0.4 0.3 0.1]/0.5;
tmp=1:1:64;  
L=6;
Ntmp=1:6;
for ntmp=1:(1*length(Ntmp)/3)
%     confounds(ntmp,:)=sin(ntmp/4*pi*tmp);
    confounds(ntmp,:)=((2/N)^0.5)*cos((2*pi/3.7)*Ntmp(ntmp)*tmp);
end
for ntmp=(1*length(Ntmp)/3)+1:(2*length(Ntmp)/3)
%     confounds(ntmp,:)=sin(ntmp/4*pi*tmp);
    confounds(ntmp,:)=((2/N)^0.5)*cos((2*pi/3.4)*Ntmp(ntmp)*tmp);
end
for ntmp=(2*length(Ntmp)/3)+1:(3*length(Ntmp)/3)
%     confounds(ntmp,:)=sin(ntmp/4*pi*tmp);
    confounds(ntmp,:)=((2/N)^0.5)*cos((2*pi/3.0)*Ntmp(ntmp)*tmp);
end
RCorr=diag(zeros(L,1)+1);
% %%%% 1st order sub-diagonal matrix%%%%%%%%%%    
% for m=1:L-1
%         RCorr(m+1,m)=0.5;
%         RCorr(m,m+1)=0.5;
%     end
% %%%%%%%%  %%%%%%%%% 

%%%%%% 2nd order sub-diagonal matrix%%%%%%%    
   for m=1:L-1
        RCorr(m+1,m)=0.5;
        RCorr(m,m+1)=0.5;
    end     
    for m=1:L-2
        RCorr(m+2,m)=0.2;
        RCorr(m,m+2)=0.2;
    end 
RCorrSqrt=sqrtm(RCorr);
confoundsR=RCorrSqrt*confounds; %correlated confounds
% SignalConf=rn*confoundsR;% final confounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% %%%% expected confounds matrix%%%%%%%%%%   
for ntmp=1:length(Ntmp)
%     confounds(ntmp,:)=sin(ntmp/4*pi*tmp);
    confounds(ntmp,:)=((2/N)^0.5)*cos((2*pi/3.7)*Ntmp(ntmp)*tmp);
end
%%%%%%%%%%%%%%%%%%%%% 
%==================================================================


%%%%%%%%%%%%%% III. This part is to load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% and change some parmeters in DCM from the loaded data.%%%%%%

SCK=load_DEM_demo_hdm_SCK();

% load SCK SCK;

DCM.Y=SCK.Y+rn*confoundsR; %%% old edition: remove _R; 20190314 edition: add _R
T=length(DCM.Y);
DCM.pU.v=SCK.pU.v{1,2};
DCM.M(2).v=SCK.pU.v{1,2};
DCM.pP.P=SCK.pP.P{1,1};
% DCM.M(1).f=SCK.M(1).f;
DCM.M(1).f = fcnchk(SCK.M(1).f,'x','v','P');
% DCM.M(1).g=SCK.M(1).g;
SCK.M(1).g='spm_gx_hdm_sck_X0';
DCM.M(1).g = fcnchk(SCK.M(1).g,'x','v','P');
DCM.M(1).ip=SCK.M(1).ip;

SCK.M(1).pE(SCK.M(1).ip)=0.2;
DCM.M(1).pE=SCK.M(1).pE;

tmp=zeros(1,7);
tmp(SCK.M(1).ip)=0.1;
DCM.M(1).pC=diag(tmp);

tmp=zeros(1,7);
tmp(SCK.M(1).ip)=0.01;
DCM.M(1).wP=diag(tmp);

tmp=zeros(1,7);
tmp(SCK.M(1).ip)=6e-9;
DCM.M(1).pP=diag(tmp);

DCM.cf.cf=confounds;
tmp=zeros(1,6)+6e-9;
DCM.cf.pP=diag(tmp);
tmp=zeros(1,6)+0.01;
DCM.cf.wP=diag(tmp);


%%%%%%%%%%%%%% IV. This part is to do model inversion %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% and the model is for confounds.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================

FULL = spm_SCKS_sDCM_X0(DCM);
%--------------------------------------------------------------------------


%%%%%%%%%%%%%% V. This part is Martin's demo %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% where results are plot.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================
%  Plot results for single time course estimates
%==================================================================
% Display State results:
f5 = spm_figure('Create','Graphics','CKF estimates');
set(f5,'RendererMode','auto','Renderer','painter');
clf(f5);
for p = 1:3
    subplot(3,1,p),
    hax = gca;
    si    = spm_invNcdf(1 - 0.05);

    if p == 1,
        xxfig = FULL.qU.r{1}(1:n,1:end);
        Sfig  = zeros(n,size(xxfig,2));
        tit   = 'BOLD prediction';
    elseif p == 2
        xxfig = FULL.qU.x{1}(1:n,1:end);
        Sfig  = FULL.qU.S{1}(1:n,1:end);
        tit   = 'Neural Estimate';
    else
        xxfig = FULL.qU.x{1}(n+1:end,1:end);
        Sfig  = FULL.qU.S{1}(n+1:end,1:end);
        tit   = 'States Estimate';
    end

    s = abs(Sfig);

    % conditional covariances
    %------------------------------------------------------------------

    j       = [1:size(xxfig,1)];
    ss      = si*s(j,:);
    s(j,:)  = [];
    [ill indss] = sort(mean(ss,2),'descend');
    pf = plot(1:T,xxfig,'linewidth',1.5);
    set(hax,'xlim',[1,T],'nextplot','add')
    for ic = 1:size(xxfig,1)
        col0 = get(pf(indss(ic)),'color');
        col = (ones(1,3)-col0)*0.65 + col0;
        fill([(1:T) fliplr(1:T)],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
            'r',...
            'FaceColor',col,...
            'EdgeColor',col);

        hold on
        COL{ic} = col0;
    end
    for ic = 1:size(xxfig,1)
        plot(xxfig(indss(ic),:),'color',COL{ic},'linewidth',0.75);
    end
    try
     if p == 1   
        plot(FULL.Y','Color',[0 0 0],'linewidth',1);
        try
          plot(SIM.pU.r{1}(1:n,:)' - SIM.pU.z{1}(1:n,:)' ,'Color','r','linewidth',1);    
        catch
        end
     elseif p == 2
        plot(SIM.pU.x{1}(1:n,:)','Color',[0 0 0],'linewidth',1);   
     else
        plot(SIM.pU.x{1}(n+1:end,:)','Color',[0 0 0],'linewidth',1);    
     end
    catch
    end
    title(tit);
    grid(hax,'on')
    axis(hax,'tight')
    set(hax,'box','on','Layer','top');
    set(hax,'tickdir','out')
   
end
    