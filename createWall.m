function [AeffW,wallR,LmxMin,LmxMax,LmyMin,LmyMax,LmzMin,LmzMax,d1W,d2W,d3W,cntW] = createWall(fLNum,dsA,dsB,ya,za,yb,zb,boxSize,xa,xb,AeffA,AeffB)

    wallNum = linspace(1,round(1.5.*fLNum),round(1.5.*fLNum))';
    [rwy,rwz] = meshgrid(wallNum,wallNum);
    rwy = rwy.*max(dsA,dsB);
    rwz = rwz.*max(dsA,dsB);
    rwy = rwy + mean(mean(ya))-mean(mean(rwy));
    rwz = rwz + mean(mean(za))-mean(mean(rwz));
    dsW = rwy(1,2)-rwy(1,1);
    AeffW = 3.5.*dsW;
    rwx = zeros(round(1.5.*fLNum),round(1.5.*fLNum))+boxSize./2 - 2.*AeffA -AeffW;
    
    
    rlstX = reshape(rwx,[size(rwx,1).*size(rwx,2),1]);
    rlstY = reshape(rwy,[size(rwy,1).*size(rwy,2),1]);
    rlstZ = reshape(rwz,[size(rwz,1).*size(rwz,2),1]);
    
    wallR = [rlstX,rlstY,rlstZ];
    cntW = size(wallR,1);
    
    d1W = [1,0,0].*ones(cntW,1);
    d2W = [0,1,0].*ones(cntW,1);
    d3W = [0,0,1].*ones(cntW,1);
    
    LmxMin = min(min(xa),min(xb))-2.0.*max(AeffA,AeffB);
    LmxMax = max(max(xa),max(xb))+10.*max(AeffA,AeffB);
    
    LmyMin = ya(1)-10.*max(AeffA,AeffB);
    LmyMax = ya(end)+10.*max(AeffA,AeffB);
    
    LmzMin = min(za(1),zb(1))-3.*max(AeffA,AeffB);
    LmzMax = min(za(end),zb(end))+3.*max(AeffA,AeffB);
end
