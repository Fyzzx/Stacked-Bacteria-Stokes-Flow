function f=vec_spread3D(F,r3D,hFL,fLNum,GridNum)
% spread F to grid

c = 1/(hFL*hFL*hFL);
f=zeros(fLNum,fLNum,fLNum,3);

s = r3D/hFL; % Get body position relative to grid
i = floor(s);
rdiff = s-i;
w = vec_phi1(rdiff(:,1)).*vec_phi2(rdiff(:,2)).*vec_phi3(rdiff(:,3));%Evaluate delta function
w = permute(w, [1,3,4,2]);

for k=1:GridNum
  i1=mod((i(k,1)-1):(i(k,1)+2),fLNum)+1; %Find affected cells
  i2=mod((i(k,2)-1):(i(k,2)+2),fLNum)+1;
  i3=mod((i(k,3)-1):(i(k,3)+2),fLNum)+1;
  ww = w(:,:,:,k);
  f(i1,i2,i3,1)=f(i1,i2,i3,1)+(c*F(k,1))*ww; %Spread force to fluid
  f(i1,i2,i3,2)=f(i1,i2,i3,2)+(c*F(k,2))*ww;
  f(i1,i2,i3,3)=f(i1,i2,i3,3)+(c*F(k,3))*ww;
end
end