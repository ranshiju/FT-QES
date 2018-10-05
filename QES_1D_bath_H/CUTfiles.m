Pwd='.\'; Kind='mat';
L=12;
if(~exist([Pwd,'DataL',num2str(L)],'dir'))
    mkdir([Pwd,'DataL',num2str(L)])
end

DataEXP=['DMRGv5J(0,1)-(0,1)_h(*,0)-(*,0)Jb1_dc(2,10)_L(4,',num2str(L),')_tau1e-08'];
[Name,~,~,~,~]=SearchFiles(Pwd,[DataEXP,'.',Kind],0);
load([Pwd,Name{1}],'paraDMRG')

DataEXP=['DMRGv5J(0,1)-(0,1)_h(',num2str(paraDMRG.hx),',0)-(*,0)Jb1_dc(2,10)_L(4,',num2str(L),')_tau1e-08'];
[Name,~,~,~,Num]=SearchFiles(Pwd,[DataEXP,'.',Kind],0);
hbath=zeros(1,Num);
for n=1:Num
    load([Pwd,Name{n}],'para')
    hbath(n)=para.hx;
end

Ni=numel(hbath);

for ni=1:Ni
    DataEXP=['DMRGv5J(0,1)-(0,1)_h(*,0)-(',num2str(hbath(ni)),',0)Jb1_dc(2,10)_L(4,',num2str(L),')_tau1e-08'];
    [Name,~,~,~,Num]=SearchFiles(Pwd,[DataEXP,'.',Kind],0);
    for n=1:Num
        movefile([Pwd,Name{n}],[Pwd,'DataL',num2str(L),'\']);
    end
end

fprintf([num2str(Ni),'x',num2str(Num),'files have been moved to "DataL',num2str(L),'". \n']);