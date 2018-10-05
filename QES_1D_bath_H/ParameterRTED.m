function para=ParameterRTED
para.L=8; % Total length, bath included, must be EVEN

para.Jbath=1;
para.Jxy=0;
para.Jz=2;
para.hx=10;
para.hz=0;
para.time=[0.01,0.05,0.1:0.2:2,2.5:0.5:10,11:40,42:2:60]; % Do NOT contain time=0 here;

% For calculating the bath Hamiltonian after quenching
para.tau0=0.1; 
para.dtau=0.1;
para.taut=8;
para.Eps=1e-10;
para.timeGS0=10;
para.timeGS=10000;
para.timeGS2=1;

para.tau=1e-4;
para.tau1=para.tau0*para.dtau^(para.taut-1);
% para.const=0;
% [para]=Check_para(para);
end