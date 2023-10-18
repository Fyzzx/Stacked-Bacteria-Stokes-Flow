function [ra,rb,xa,ya,za,xb,yb,zb,d1a,d1b,d2a,d2b,d3a,d3b,rOriga,rOrigb] = initLocs(aA,aB,sA,sB,phiA,phiB,GNA,GNB,boxSize,AeffA,AeffB,yshift,zshift)
    ra(:,:) = [aA  ,sA'.*sin(phiA), sA'.*cos(phiA)]; % x dimension is height off wall
    rb(:,:) = [aA(1)+2*aB,sB'.*sin(phiB), sB'.*cos(phiB)];
    
    xa = ra(:,1);
    ya = ra(:,2);
    za = ra(:,3);
    
    xb = rb(:,1);
    yb = rb(:,2);
    zb = rb(:,3);
    
    d1a = [1,0,0].*ones(GNA,1);
    d1b = [1,0,0].*ones(GNB,1);
    d2a = [0,cos(phiA),-sin(phiA)].*ones(GNA,1);
    d2b = [0,cos(phiB),-sin(phiB)].*ones(GNB,1);
    d3a = [0,sin(phiA),cos(phiA)].*ones(GNA,1);
    d3b = [0,sin(phiB),cos(phiB)].*ones(GNB,1);
    
    xa = xa + boxSize./2 - 2*AeffA;
    ya = ya - mean(ya);
    ya = ya + boxSize./2;
    za = za - mean(za);
    za = za + boxSize./2 - mean(za);
    
    xb = xb + boxSize./2 - AeffA - AeffB;
    if abs(phiA-phiB) < pi/16
        xb = xb + ones(GNB,1).*(sqrt((2*AeffB).^2-(yshift).^2)-2.*AeffB);
    end
    yb = yb -mean(yb);
    yb = yb + boxSize./2  + yshift;
    zb = zb - mean(zb);
    zb = zb + boxSize./2  + zshift;
    
    rOriga = [xa,ya,za];
    rOrigb = [xb,yb,zb]; % this sets up the bacteria to not move
end
