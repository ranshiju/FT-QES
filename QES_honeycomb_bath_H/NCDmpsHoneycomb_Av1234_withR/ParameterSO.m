function para=ParameterSO
% 物理参数
para.Jh=1;
para.Jk=0;
para.hx=0;
para.hz=0;
para.tau0=0.1;
para.dtau=0.1;
para.taut=6;

% 计算参数
para.dc=8;
para.tau=[1e-1,1e-2,1e-3,1e-4,1e-5];
para.Eps=1e-8;
para.time=1000;
para.time0=20;
para.dtime=5; % The interval to check the convergence

para.dt_print=20; % The interval to print the results
para.IsHermi=1;
para.IsSymme=1; % If v1=v3; v2=v4
para.IsReal=1;
end