chi=6;
d=2;

T1=randn(d,chi,chi,chi);
T2=randn(d,chi,chi,chi);
LM1=sort(rand(chi,1),'descend');
LM2=sort(rand(chi,1),'descend');
LM3=sort(rand(chi,1),'descend');

[T,LM]=SuperOrthogonalizationHoneycomPEPS({T1,T2},{LM1,LM2,LM3},chi);
Ort=TestSuperOrthogonality(T,LM);
fprintf('The deviation from SO form is %g \n',Ort)