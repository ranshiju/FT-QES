function [v,S]=v_unitary(v,D,dc,Way)
if(Way==1)
    v.L=reshape(v.L,[dc,D*dc]);
    v.R=reshape(v.R,[dc,D*dc]);

    [L,S,R]=svd(v.L.'*v.R); S=diag(S);

    v.L=reshape( L(:,1:dc)', [dc*D*dc,1] );
    v.R=reshape( R(:,1:dc).',[dc*D*dc,1] );
elseif(Way==2)
    dh=round(sqrt(D));
    v.L=reshape( permute( reshape(v.L,[dc,dh,dh,dc]),[1,2,4,3] ),[dc*dh,dc*dh] );
    v.R=reshape( permute( reshape(v.R,[dc,dh,dh,dc]),[1,2,4,3] ),[dc*dh,dc*dh] );
    
%     [v.L,v.R]=TruncationFromRM(v.L,v.R,ones(dc*dh,1),0);
%     reshape(v.L(:,1:dc),[])
end
end