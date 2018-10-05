function para=Parameter
para.ISDMRG=0;
para.ISED=0;
para.ISFTED=0;
para.ISRTED=0;
para.ISLTRG=1;

para.Jxy=1;
para.Jz=0; % For real-time evolution, set Jz=2 to compare with the analytical solution
para.hx=0;
para.hx2=para.hx; % The x-field on the second and last second sites
para.hz=0;
para.L=4; % Even number, bulk size

para.dc=2;
para.tau=1e-1;
para.taut=8;
para.dtau=0.1;
para.Nstate=2;
para.ISeig=1;
para.Otrotter=1; % Trotter order

para.time=10000;
para.time0=10;
para.time_canon=10000;
para.time2=1;
para.Eps=1e-10;

para.ISsparse=1;
para.tau1=para.tau*para.dtau^(para.taut-1);
[para]=Check_para(para);
end