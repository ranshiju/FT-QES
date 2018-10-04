function para=Parameter
para.ISDMRG=0; % DMRG for ground state
para.ISFTED=0; % ED for finite-T
para.ISLTRG=1; % MPO-LTRG for finite-T

% Physical parameters
para.Jh=1;
para.Jk=0;
para.hx=0;
para.hz=0;

% Computational parameters for physical-bath interactions
para.tau0=0.1;
para.dtau=0.1;
para.taut=6;

para.dc=2; % Dimension of the bath sites
para.time=3000;
para.dtime=10;
para.time0=20;
para.dt_print=50;
para.Eps=1e-9; % Threshold to break the loop

% Computational parameters in the second step (DMRG for QES)
para.ISsvd=0; % svd or QR
para.Conv=1e-9;
para.NpMPS=[2,5,6,7,10,11,12,15]; % number of physical sites in MPS
para.NbMPS=[1,3,4,8,9,13,14,16]; % number of bath sites in MPS

% Parameter for DMRG
para.chi=8; % Dimension of the DMRG MPS
para.timeDMRG=100;
para.tauDMRG=para.tau0*(para.dtau^(para.taut-1));
para.Nmps=numel(para.NpMPS)+numel(para.NbMPS);
para.N0=2; % How many physical sites on one tensor

para.IsHermi=1;
para.IsSymme=1; % If v1=v3; v2=v4
para.IsReal=1;

% eigs parameters
para.p=6;
para.maxit=10000;
para.tol=1e-12;
para.isreal=1;
end