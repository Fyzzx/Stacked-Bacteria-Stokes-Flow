function [bDelta] = init_bDelta_3D(fLNum,hFL,dt,rho)
%init_b.m
%This script initializes the array  b 
%that is used in the fluid solver to extract pressure output

bDelta=zeros(fLNum,fLNum,fLNum,3);

for m1=0:(fLNum-1)
  for m2=0:(fLNum-1)
    for m3=0:(fLNum-1)
        if~(((m1==0)||(m1==fLNum/2))&&((m2==0)||(m2==fLNum/2))&&((m3==0)||(m3==fLNum/2)))
          t=2*pi/fLNum*[m1;m2;m3];
          s = sin(t);
          c = - (1i*rho*hFL/dt)/(s(1)^2+s(2)^2+s(3)^2);
          bDelta(m1+1,m2+1,m3+1,1) = c*s(1);
          bDelta(m1+1,m2+1,m3+1,2) = c*s(2);
          bDelta(m1+1,m2+1,m3+1,3) = c*s(3);
        end
    end
  end
end
end