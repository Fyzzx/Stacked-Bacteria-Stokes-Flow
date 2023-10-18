function [uuu]=fluid3D_IB(uFL,fUpFL,aDelta,dt,eta,ipIB,imIB,hFL,rho)
% function [uuu,uu,ppp]=fluid3D_IB(uFL,fUpFL,aDelta,bDelta,dt,eta,ipIB,imIB,hFL)

%uuu outputs the final velocity, where uu is the output in between. 
    %See notes on https://www.math.nyu.edu/faculty/peskin/ib_lecture_notes/index.html 
    %ppp outputs the final pressure.
%     w=uFL-(dt/2)*skew(uFL)+(dt/(2*rho))*fUpFL; %includes (u dot Del)u
%     w = uFL+(dt/(2*rho))*fUpFL;
% 
%     w=fft(w,[],1);
%     w=fft(w,[],2);
%     w=fft(w,[],3);
%     
%     uu(:,:,:,1)=aDelta(:,:,:,1,1).*w(:,:,:,1)+aDelta(:,:,:,1,2).*w(:,:,:,2)+aDelta(:,:,:,1,3).*w(:,:,:,3);
%     uu(:,:,:,2)=aDelta(:,:,:,2,1).*w(:,:,:,1)+aDelta(:,:,:,2,2).*w(:,:,:,2)+aDelta(:,:,:,2,3).*w(:,:,:,3);
%     uu(:,:,:,3)=aDelta(:,:,:,3,1).*w(:,:,:,1)+aDelta(:,:,:,3,2).*w(:,:,:,2)+aDelta(:,:,:,3,3).*w(:,:,:,3);
%     
%     
%     uu=ifft(uu,[],3);
%     uu=ifft(uu,[],2);
%     uu=real(ifft(uu,[],1));
    
%     w=uFL-dt*skew(uu)+(dt/rho)*fUpFL+(dt/2)*(mu/rho)*laplacian3D(uFL,imIB,ipIB,hFL); %includes (u dot Del)u
    w=uFL+(dt/rho)*fUpFL+(dt/2)*(eta/rho)*laplacian3D(uFL,imIB,ipIB,hFL);

    w=fft(w,[],1);
    w=fft(w,[],2);
    w=fft(w,[],3);
    
    
    uuu(:,:,:,1)=aDelta(:,:,:,1,1).*w(:,:,:,1)+aDelta(:,:,:,1,2).*w(:,:,:,2)+aDelta(:,:,:,1,3).*w(:,:,:,3);
    uuu(:,:,:,2)=aDelta(:,:,:,2,1).*w(:,:,:,1)+aDelta(:,:,:,2,2).*w(:,:,:,2)+aDelta(:,:,:,2,3).*w(:,:,:,3);
    uuu(:,:,:,3)=aDelta(:,:,:,3,1).*w(:,:,:,1)+aDelta(:,:,:,3,2).*w(:,:,:,2)+aDelta(:,:,:,3,3).*w(:,:,:,3);
    
%     ppp=bDelta(:,:,:,1).*w(:,:,:,1)+bDelta(:,:,:,2).*w(:,:,:,2)+bDelta(:,:,:,3).*w(:,:,:,3);
%     ppp=ifft(ppp,[],3);
%     ppp=ifft(ppp,[],2);
%     ppp=real(ifft(ppp,[],1));
    %ppp is the pressure output
    uuu=ifft(uuu,[],3);
    uuu=ifft(uuu,[],2);
    uuu=real(ifft(uuu,[],1));
end