function para=ParameterRC
% Parameters of the bulk for remote control
para.L=15; % The true lenth is L*2 + 2 bath sites
para.Jxy=0;
para.Jz=1;
para.hx=0.5;
para.hz=0;

% Bath coupling strength
para.Jbath=1; % Hbath=Jbath*Hbath

% DMRG parameters
para.tauDMRG=1e-4;
para.chi=16;
para.timeDMRG=100;
para.isreal=1;
para.Jeps=1e-15;

para.N0=2; % Numer of spins per tensor
para.IsObBath=1; % Allowed only with bath dimension 2
end