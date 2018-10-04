function MainTree(varargin)
%% Preperation
fprintf('Start simulation of the physical-bath interactions \n');
d=2;
Ni=numel(varargin);
if(Ni==0)
    paraTree=ParameterTree;
    paraTree.PATHtree=[pwd,filesep];
    paraTree.EXP=['J_(',num2str(paraTree.Jh),',',num2str(paraTree.Jk),')_h_(',num2str(paraTree.hx),',',num2str(paraTree.hz),')_dc',num2str(paraTree.dc),...
        '_tau',num2str(paraTree.tau0*(paraTree.dtau^(paraTree.taut-1))),'_ISSymme',num2str(paraTree.IsSymme)];
elseif(Ni==1)
    paraTree=varargin{1};
end

%% Initalize state and define Hamiltonian
[Hx,Hy,Hz]=HamiltonianKitaHesen(paraTree.Jh,paraTree.Jk,paraTree.hx,paraTree.hz);
A=InitialPEPS(paraTree.dc,d);
[Vx]=InitialEnv(paraTree.dc,d^2); Vy=Vx;

paraTree.tau=paraTree.tau0*(paraTree.dtau.^(0:paraTree.taut-1));
Ntau=paraTree.taut;
Ez=zeros(paraTree.time,1); Ex=Ez; Ey=Ez; Etot=Ez;
OBsimp.Ez=zeros(Ntau,1); OBsimp.Ex=OBsimp.Ez; OBsimp.Ey=OBsimp.Ez;
OBsimp.Ef=zeros(Ntau,1);
Size=[d,d,paraTree.dc,paraTree.dc,paraTree.dc,paraTree.dc];

%%Iteration over tau
for nt=1:Ntau
    fprintf('Simulation with tau = %g  \n',paraTree.tau(nt));
    %% Prepare the evolution gates
    [ULx,URx]=ExpSvdU(Hx,paraTree.tau(nt));
    [ULy,URy]=ExpSvdU(Hy,paraTree.tau(nt));
    Ub=expm(-paraTree.tau(nt)*Hz);
    tob=0;
    for t=1:paraTree.time
        %% Iteration with each fixed tau
        IsOb=(rem(t,paraTree.dtime)==0);
        opts.v0=A(:);
        [A,~]=eigs(@(A) fhandle_updateA(A,Ub,{ULx,ULy,URx,URy},{Vx,Vy},Size,paraTree.IsHermi),d*d*paraTree.dc^4,1,'lm',opts);
        A=reshape(A,Size);
        if(paraTree.IsReal)
            A=real(A);
        end
        if(IsOb)
            tob=tob+1;
            [Ez(tob)]=Observation(A,Hz);
        end
        
        [Aort]=Orthogonalize_A(A,5);
        if(IsOb)
            [Ex(tob)]=ObservationXY(A,Aort,Hx,1);
        end
        Vx=UpdateVxy(Aort,ULx,ULy,URx,URy,Ub,Vx,Vy,5,paraTree.IsHermi);
        
        [Aort]=Orthogonalize_A(A,6);
        if(IsOb)
            [Ey(tob)]=ObservationXY(A,Aort,Hy,2);
            Etot(tob)=(Ex(tob)+Ey(tob)+Ez(tob))/2;
        end
        Vy=UpdateVxy(Aort,ULx,ULy,URx,URy,Ub,Vx,Vy,6,paraTree.IsHermi);
        
        IsBreak=(t>paraTree.time0 && tob>1.1 && abs(Etot(tob)-Etot(tob-1))<paraTree.Eps);
        if(IsBreak)
            fprintf('Converged at the %g-th iteraction, Etot = %g \n',t,Etot(tob))
            break;
        elseif(rem(t,paraTree.dt_print)==0)
            fprintf('At the %g-th iteration, Etot = %g \n',t,Etot(tob))
        end
    end
    OBsimp.Ex(nt)=Ex(tob);
    OBsimp.Ey(nt)=Ey(tob);
    OBsimp.Ez(nt)=Ez(tob);
    OBsimp.Ef(nt)=Etot(tob);
end

save([paraTree.PATHtree,paraTree.EXP,'.mat'],'paraTree','A','Vx','Vy','ULx','URx','ULy','URy','Ub','Hx','Hy','Hz','OBsimp');
end

function A=InitialPEPS(chi,d)
A=randn(d,d,chi,chi,chi,chi);
A=A+permute(A,[2,1,3,4,5,6]);
A=A+permute(A,[1,2,3,6,5,4]);
A=A+permute(A,[1,2,5,4,3,6]);
end

function [V]=InitialEnv(chi,D)
V=randn(chi,D,chi);
V=(V+permute(V,[3,2,1]))/2;
end

function [UL,UR]=ExpSvdU(H,tau)
d=2;
% U=reshape( permute( reshape(expm(-tau*H),[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
U=reshape( permute( reshape(eye(size(H))-tau*H,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );

[UL,LM,UR]=svd(U);
UL=UL*sqrt(LM); 
UR=conj(UR)*sqrt(LM);
end

function A1=fhandle_updateA(A,Ub,Uh,V,Size,IsHermi)
A=reshape(A,Size);

A1=TwoBodyEvolution(Ub,A,[1,2]);
U2=LocalH_UV(Uh{3},V{1},IsHermi);
A1=TwoBodyEvolution(U2,A1,[1,3]);
U2=LocalH_UV(Uh{4},V{2},IsHermi);
A1=TwoBodyEvolution(U2,A1,[1,4]);
U2=LocalH_UV(Uh{1},V{1},IsHermi);
A1=TwoBodyEvolution(U2,A1,[2,5]);
U2=LocalH_UV(Uh{2},V{2},IsHermi);
A1=TwoBodyEvolution(U2,A1,[2,6]);

A1=A1(:);
end

function A1=TwoBodyEvolution(H2,A,p)
S=size(A); Nb=numel(S);

A1=permute(A,[p(1),p(2),1:p(1)-1,p(1)+1:p(2)-1,p(2)+1:Nb]);
A1=H2*reshape(A1,[S(p(1))*S(p(2)),prod(S([1:p(1)-1,p(1)+1:p(2)-1,p(2)+1:Nb]))]);
A1=permute( reshape(A1,S([p(1),p(2),1:p(1)-1,p(1)+1:p(2)-1,p(2)+1:Nb])),[3:p(1)+1,1,p(1)+2:p(2),2,p(2)+1:Nb] );
end

function H2=LocalH_UV(U,V,IsHermi)
[dd,D]=size(U); chi=size(V,1);
d=round(sqrt(dd));

V=permute(V,[2,1,3]);
H2=reshape( permute( reshape( U*reshape(V,[D,chi*chi]),[d,d,chi,chi] ),[1,3,2,4] ),[d*chi,d*chi] );
if(IsHermi)
    H2=(H2+H2')/2;
end
end

function A=SymmetrizeA(A) %#ok<DEFNU>
A=(A+permute(A,[1,2,5,4,3,6]))/2;
A=(A+permute(A,[1,2,3,6,5,4]))/2;
end

function [Ap]=Orthogonalize_A(A,p)
[d,~,chi,~,~,~]=size(A);

Ap=reshape( permute(A,[1:p-1,p+1:6,p]),[d*d*chi^3,chi] );
[Ap,~]=qr(Ap,0);
Ap=permute( reshape(Ap,[d,d,chi,chi,chi,chi]),[1:p-1,6,p:5] );
end

function Vp=UpdateVxy(Ap,ULx,ULy,URx,URy,Ub,Vx,Vy,p,IsHermi)
[d,~,chi,~,~,~]=size(Ap);
D=d*d;
A1=TwoBodyEvolution(Ub,Ap,[1,2]);
U2=LocalH_UV(URx,Vx,IsHermi);
A1=TwoBodyEvolution(U2,A1,[1,3]);
U2=LocalH_UV(URy,Vy,IsHermi);
A1=TwoBodyEvolution(U2,A1,[1,4]);

A1=reshape( permute(A1,[1,3,4,2,5,6]),[d*chi*chi,d*chi*chi] );
A1=reshape( permute(conj(Ap),[2,5,6,1,3,4]),[d*chi*chi,d*chi*chi] )*A1;
A1=reshape(A1,[d,chi,chi,d,chi,chi]);
A1=reshape( permute(A1,[2,3,5,6,1,4]),[chi^4,d*d] );
Vp=reshape( permute( reshape( reshape( permute( reshape(ULy,[d,d,D]),[1,3,2] ),[d*D,d] ) * reshape(ULx,[d,d*D]),[d,D,d,D] ),[1,3,4,2] ),[d*d,D*D] );
A1=reshape( permute( reshape( A1 * Vp,[chi,chi,chi,chi,D,D] ),[1,5,3,2,6,4] ),[chi*D*chi,chi*D*chi] );

if(p==5) % Update Vx
    Vp=reshape( A1*reshape(Vy,[chi*D*chi,1]),[chi,D,chi] );
elseif(p==6) % Update Vy
    Vp=reshape( reshape(Vx,[1,chi*D*chi])*A1,[chi,D,chi] );
end
Vp=Vp/norm(Vp(:));
end

function rho2=RDM_TwoBody(A)
[d,~,d1,d2,d3,d4]=size(A);
A=reshape(A,[d^2,d1*d2*d3*d4]);
rho2=A*A';
rho2=rho2/trace(rho2);
end

function [E]=Observation(A,H)
    Rho2=RDM_TwoBody(A);
    E=trace(H*Rho2);
end

function [Exy]=ObservationXY(A,Aort,H,n)
% Always observe on the 1st and 2nd virtual bonds
% n=1 or 2
[d,~,d1,d2,d3,d4]=size(A);
if(n==1)
    A=reshape( permute(A,[2,4,5,6,1,3]),[d*d2*d3*d4,d*d1] );
    A=reshape( permute( reshape( A'*A,[d,d1,d,d1] ),[2,4,1,3] ),[d1*d1,d*d] );
    
    Aort=reshape( permute(Aort,[1,3,4,6,2,5]),[d*d1*d2*d4,d*d3] );
    Aort=reshape( permute( reshape( Aort'*Aort,[d,d3,d,d3] ),[1,3,2,4] ),[d*d,d3*d3] );
elseif(n==2)
    A=reshape( permute(A,[2,3,5,6,1,4]),[d*d1*d3*d4,d*d2] );
    A=reshape( permute( reshape( A'*A,[d,d2,d,d2] ),[2,4,1,3] ),[d2*d2,d*d] );
    
    Aort=reshape( permute(Aort,[1,3,4,5,2,6]),[d*d1*d2*d3,d*d4] );
    Aort=reshape( permute( reshape( Aort'*Aort,[d,d4,d,d4] ),[1,3,2,4] ),[d*d,d4*d4] );
end
rho=reshape( permute( reshape( Aort*A,[d,d,d,d] ),[1,3,2,4] ),[d*d,d*d] );
rho=rho/trace(rho);
Exy=trace(rho*H);
end

