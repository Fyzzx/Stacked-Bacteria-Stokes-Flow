function [d1new,d2new,d3new] = updateDHat(GN,dt,omegaMag,d1,d2,d3,omegaSc)
        for mm = 1:GN % cell A
        eTerm = omegaSc(mm,:)';
        eTermT = transpose(eTerm);
        term2eeT = eTerm*eTermT;

        rotMat1(:,:) =  [cos(omegaMag(mm).*dt) ,0 ,0;...
                            0, cos(omegaMag(mm).*dt) ,0;...
                            0, 0, cos(omegaMag(mm).*dt)] ;

         
        rotMat2(:,:) = (1-cos(omegaMag(mm).*dt)).*(term2eeT);
       
    % 3rd term is the sin(theta)(e(cross))
        rotMat = rotMat1 + rotMat2;
  
        for bb = 1:3
            dlist = [d1(mm,:);d2(mm,:);d3(mm,:)];
            dnow = dlist(bb,:);
            
        [eXOldDx,eXOldDy,eXOldDz] = crossProduct(eTerm',dnow);
        eXOldD = [eXOldDx;eXOldDy;eXOldDz];
        rotMat3withD = sin(omegaMag(mm).*dt).*eXOldD;
    
        newnewDvect = rotMat(:,:)*dnow' + rotMat3withD;
        dnew(:,bb) = newnewDvect;
       
        end
           
    %% Update the bacteria bases
        d1new(mm,:) = dnew(:,1)';
        d2new(mm,:) = dnew(:,2)';
        d3new(mm,:) = dnew(:,3)';
        
        end 
end 

function [xterm,yterm,zterm] = crossProduct(term1,term2)
    xterm = term1(:,2).*term2(:,3) - term1(:,3).*term2(:,2);
    yterm = term1(:,3).*term2(:,1) - term1(:,1).*term2(:,3);
    zterm = term1(:,1).*term2(:,2) - term1(:,2).*term2(:,1);
end