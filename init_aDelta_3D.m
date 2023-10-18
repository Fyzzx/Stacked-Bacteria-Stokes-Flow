function [aDelta] = init_aDelta_3D(flNum,dt,eta,hFL,rho)
%init_a.m
%This script initializes the array  a  
%that is used in the fluid solver
mu = eta;
aDelta=zeros(flNum,flNum,flNum,3,3);
for m1=0:(flNum-1)
  for m2=0:(flNum-1)
      for m3=0:(flNum-1)
    aDelta(m1+1,m2+1,m3+1,1,1)=1;
    aDelta(m1+1,m2+1,m3+1,2,2)=1;
    aDelta(m1+1,m2+1,m3+1,3,3)=1;
      end
  end
end

for m1=0:(flNum-1)
  for m2=0:(flNum-1)
    for m3=0:(flNum-1)
        if~(((m1==0)|(m1==flNum/2))&((m2==0)|(m2==flNum/2))&((m3==0)|(m3==flNum/2)))
          t=(2*pi/flNum)*[m1;m2;m3];
          s=sin(t);
          ss=(s*s')/(s'*s);
    %     a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)-(s*s')/(s'*s);
          aDelta(m1+1,m2+1,m3+1,1,1)=aDelta(m1+1,m2+1,m3+1,1,1)-ss(1,1);
          aDelta(m1+1,m2+1,m3+1,1,2)=aDelta(m1+1,m2+1,m3+1,1,2)-ss(1,2);
          aDelta(m1+1,m2+1,m3+1,1,3)=aDelta(m1+1,m2+1,m3+1,1,3)-ss(1,3); 
          
          aDelta(m1+1,m2+1,m3+1,2,1)=aDelta(m1+1,m2+1,m3+1,2,1)-ss(2,1);
          aDelta(m1+1,m2+1,m3+1,2,2)=aDelta(m1+1,m2+1,m3+1,2,2)-ss(2,2);
          aDelta(m1+1,m2+1,m3+1,2,3)=aDelta(m1+1,m2+1,m3+1,2,3)-ss(2,3); 
          
          aDelta(m1+1,m2+1,m3+1,3,1)=aDelta(m1+1,m2+1,m3+1,3,1)-ss(3,1);
          aDelta(m1+1,m2+1,m3+1,3,2)=aDelta(m1+1,m2+1,m3+1,3,2)-ss(3,2);
          aDelta(m1+1,m2+1,m3+1,3,3)=aDelta(m1+1,m2+1,m3+1,3,3)-ss(3,3); 
    
        end
    end
  end
end

for m1=0:(flNum-1)
  for m2=0:(flNum-1)
      for m3=0:(flNum-1)
    t=(pi/flNum)*[m1;m2;m3];
    s=sin(t);
    aDelta(m1+1,m2+1,m3+1,:,:)=aDelta(m1+1,m2+1,m3+1,:,:)...
                    /(1+(dt/2)*(eta/rho)*(4/(hFL*hFL))*(s'*s));
  
      end
  end
end
end