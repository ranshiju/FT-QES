function MainBH(varargin)
tic
fprintf('Program begins with parameters: \n')
Ni=numel(varargin);
if(Ni==0)
    para=Parameter;
elseif(Ni==1)
    if ischar(varargin{1})
        load(varargin{1});
    else
        if numel(varargin{1})==1
            para = varargin{1}{1};
        elseif numel(varargin{1})==2
            para = varargin{1}{1};
        end
    end
else
    error('Bad input of Main.m');
end
para.dc1=para.dc;
para.Way=1;

if(para.ISeig)
    para.EXP=['Hbath_J(',num2str(para.Jxy),',',num2str(para.Jz),')_h(',num2str(para.hx),',',num2str(para.hz),')_dc',num2str(para.dc),...
        '_L',num2str(para.L),'_tau',num2str(para.tau1),'.mat'];
else
    para.EXP=['Hbath_Jx',num2str(para.Jx),'_hz',num2str(para.hz),'_dc',num2str(para.dc),'_L',num2str(para.L),'_tau',num2str(para.tau),'S','.mat'];
end

if(~exist(para.EXP,'file'))
para %#ok<NOPRT>
fprintf('Data of the entanglement bath do not exist. Launch AOP1D simulation. \n')

opt.isreal=1;

Ob.E=cell(para.taut,1);
Ob.Etot=Ob.E; Ob.EL=Ob.E; Ob.E=Ob.E;

Ob.Ent=zeros(para.taut,2);
Ob.Xi=Ob.Ent;
Ob.EntvL=zeros(para.taut,1);
Ob.EntvR=zeros(para.taut,1);
Ob.XivL=zeros(para.taut,1);
Ob.XivR=zeros(para.taut,1);
Ob.Xiv=zeros(para.taut,1);
Ob.Es=zeros(para.taut,para.Nstate);
Ob.Mx=zeros(para.taut,para.L);
Ob.Mz=Ob.Mx;

Info.Time=zeros(para.taut,1);
d=2; L=para.L; dc=para.dc;
tau=para.tau; 
Info.ErrE=zeros(para.taut,1); Ob.Emid=zeros(para.taut,1);
% lambda=ones(para.dc,1);

for n=1:para.taut
    fprintf('tau = %g \n',tau);
    Ob.E{n}=zeros(para.time,para.L-1);
    Ob.Etot{n}=zeros(para.time,1);
    Ob.EL{n}=zeros(para.time,1);
    Ob.ER{n}=zeros(para.time,1);
%     Ob.Ent{n}=zeros(para.time,2);
    Err=1; Eps=para.Eps;
    
    [G,H]=Local_gate(para.Jxy,para.Jz,para.hx,para.hz,tau,para.L,para.ISsparse,para.Otrotter);
    if(n==1)
        [A,v]=Initial_tensor(para.L,para.dc,d,para.Way);
%         A=rand(d^L*dc*dc,1);
    end
    
    for t=1:para.time
        [AL,AR,Info.lambdaL,Info.lambdaR]=Orthogonal_A_LR(A,d,dc,L);
        
        [ML]=TM_Update_v(AL,v,G.L,G.R,G.bulk,L,d,dc,para.Otrotter);
        [MR]=TM_Update_v(AR,v,G.L,G.R,G.bulk,L,d,dc,para.Otrotter);
        
        %         [Ob.E{n}(t,:),Ob.Econnect{n}(t),~,~,Info.lambda]=Observe_Energy(A,H,para.L,para.dc1);
        [Ob.E{n}(t,:),Ob.EL{n}(t),Ob.ER{n}(t),Ob.Mx(n,:),Ob.Mz(n,:)]=Observation_ALR(A,AL,AR,H,d,dc,L);
        Ob.Etot{n}(t)=sum(Ob.E{n}(t,:))/(para.L-1);
        
        for tt=1:para.time2
            v.L=ML.'*v.L; v.L=v.L/norm(v.L);
            v.R=MR*v.R; v.R=v.R/norm(v.R);
        end
        
        f=@(A) Update_A_handle(A,v,G,para);
        if(~para.ISeig)
            for tt=1:para.time2
                A=f(A); A=A/norm(A);
            end
        else
            opt.v0=real(A); [A,~]=eigs(f,d^para.L*para.dc1*para.dc1,1,'LM',opt);
        end
        
        if(t>=para.time0 && Err<Eps)
            
            Ob.E{n}=Ob.E{n}(1:t,:);
            Ob.Etot{n}=Ob.Etot{n}(1:t);
            Ob.EL{n}=Ob.EL{n}(1:t);
            Ob.ER{n}=Ob.ER{n}(1:t);
%             Ob.Ent{n}=Ob.Ent{n}(1:t,:);
            break;
        elseif(t>=2)
            Err=abs(Ob.Etot{n}(t-1)-Ob.Etot{n}(t)); 
%             [A,lambda]=Canonicalize_A(A,d,para.dc,para.L); 
%             [x1,x2]=TM_MPS_eig(A,d,para.dc,para.L);
%             keyboard
        end
    end
    Info.Time(n)=t;
    Info.Conv=Err;
    tau=tau*para.dtau;
    Ob.Emid(n)=Ob.E{n}(end,round(para.L/2));
    Ob.E{n}=Ob.E{n}';
    
    v1.L=permute( reshape(v.L,[dc,d*d,dc]),[2,1,3] ); v1.R=permute( reshape(v.R,[dc,d*d,dc]),[2,1,3] );
    [Ob.lambdavL,Ob.EntvL(n),Ob.XivL(n)]=Entanglement_infinite_MPS({v1.L},para.time_canon,10,1e-12,2);
%     [Ob.XivL(n)]=Correlation_TM_MPS({v1.L});
    [Ob.lambdavR,Ob.EntvR(n),Ob.XivR(n)]=Entanglement_infinite_MPS({v1.R},para.time_canon,10,1e-12,2);
%     [Ob.XivR(n)]=Correlation_TM_MPS({v1.R});
    
    [~,Ob.Es(n,:),Ob.Xiv(n)]=EffectiveEspectrum(v,G,H,L,d,dc,para);
    
    Ob.Ent(n,1)=Entanglement(Info.lambdaL);
    Ob.Ent(n,2)=Entanglement(Info.lambdaR);
    
    [Ob.Xi(n,:)]=Correlation_iDMRG(AL,AR,L,dc);
end
Ob.Entf=Ob.Ent(end,:);

Hpb{1}=Obtain_Hfb(v.L,G.R,para.dc,d,para.tau1);
Hpb{2}=Obtain_Hfb(G.L,v.R,d,para.dc,para.tau1);

save(para.EXP,'Hpb');
else
    fprintf('Data of the entanglement bath exist. Load directly. \n')
    load(para.EXP);
end
end