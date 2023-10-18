function [boxSize,hFL,fLNum,ipIB,imIB,ygrid,xgrid,zgrid] = initFluid(LA,LB,dsA,dsB)
    boxSize = 6.0.*max(LA,LB); % length on each side of the cube that defines the fluid
    hFL = 2*max(dsA,dsB); % sets the length of each side of the cube fluid grid
    fLNum = round(boxSize./hFL,0); % Number of fluid nodes in each cartesian dimension, make sure this is an integer by nice choice of len and GridNum
    ipIB=[(2:fLNum),1]; % Grid index shifted left
    imIB=[fLNum,(1:(fLNum-1))]; % Grid index shifted right
    meshlin = linspace(0,boxSize-hFL, fLNum);
    [ygrid,xgrid,zgrid] = meshgrid(meshlin,meshlin,meshlin); % used in quiver plot
end
