Pwd='C:\Users\sran\Dropbox\MatlabCode\NCDmpsHoneycomb_Av1234_withR\';
DataEXP='J_(1,1)_h_(0,0)_dc*_ISconj2';
% DataEXP='testJx1_hz0.5_dc*_L*_tau0.1';
Kind='mat';

% Data=LoadDataByName(Pwd,DataEXP,Kind,'dc','Ent','Xi');
[Name,~,~,~,Num]=SearchFiles(Pwd,[DataEXP,'.',Kind],0);

Ent=zeros(Num,1); x=Ent; Xi=Ent;
for n=1:Num
    load([Pwd,Name{n}],'para','Ob')
    
    tau=para.tau0*(para.dtau^(para.taut-1));
    x(n)=para.dc;
    
    Ent(n)=Entanglement(Ob.LM{1});
    Xi(n)=Ob.Xit(end)*tau;  
end
Vecs1=OrderVectors({x,Ent,Xi},'ascend');
x=Vecs1{1};
Ent=Vecs1{2};
Xi=Vecs1{3};

% subplot(1,2,1);
plot(Xi,Ent,'*red')
% subplot(1,2,2);
% plot(Xis,Ents,'oblue')