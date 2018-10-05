function para=ParameterLTRG_QES()
para.debug=0;

para.beta=[0.02,0.05, 0.1:0.1:2, 2.2:0.2:4, 4.5, 5, 6:10, 12:2:30];
para.L=40;
para.dc=400;
para.tau=0.01;

para.coeff = 1; % this satisfies E_per_site = E_per_bond * coeff (chain:1, honeycomg:3/2, square: 2)
para.ObWay=2; % 1: silayer MPO observation; 2: bilayer MPO observation
para.MPOway = 1; % 1: traditional way; 2: Tylar way (by W.Li, not allowed for QES)
end