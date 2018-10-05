function MainBH(varargin)
%% Calcualte the ground state or FT density matrix of 2D model
% This code is particularly for Heisenberg + Kitaev model
% The interactions on three bonds are inequivalent

fprintf('Start running the code: QES with entanglement bath with parameters: \n');
addpath(genpath('.'));

%% Preperation
if(numel(varargin)==0)
    para=Parameter;
elseif(numel(varargin)==1)
    load(varargin{1});
end
para.EXP=['J_(',num2str(para.Jh),',',num2str(para.Jk),')_h_(',num2str(para.hx),',',num2str(para.hz),')_dc',num2str(para.dc),'_tau',...
    num2str(para.tauDMRG),'_IsSymme',num2str(para.IsSymme)];
para %#ok<NOPRT>
para.tau1 = para.tau0 * para.dtau^(para.taut-1);

para.PATHtree=['.',filesep,'NCDmpsHoneycomb_Av1234_withR',filesep];
para.PATHDMRG=['.',filesep,'DMRGHoneycomb',filesep];
if(~exist([para.PATHtree,para.EXP,'.mat'],'file'))
    fprintf('Data of bath does no exist. Run MainTree.m  \n');
    MainTree(para);
else
    fprintf('Data of bath exists. Load directly.  \n');
end
load([para.PATHtree,para.EXP,'.mat'],'Vx','Vy','ULx','URx','ULy','URy')

%% Calculate the bath-physical interactions
fprintf('Obtain the interactions between the physical and bath sites.  \n');
[~,~, ~, Hpb{1}]=Boundary_Interactions_site(Vx,URx,para.ISsvd, para.tau1); 
[~,~, ~, Hpb{2}]=Boundary_Interactions_site(Vy,URy,para.ISsvd, para.tau1);
[~,~, ~, Hpb{3}]=Boundary_Interactions_site(Vx,ULx,para.ISsvd, para.tau1);
[~,~, ~, Hpb{4}]=Boundary_Interactions_site(Vy,ULy,para.ISsvd, para.tau1);

save(['Hbath_', para.EXP], 'Hpb')
end