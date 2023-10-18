function [w] = laplacian3D(uFL,imIB,ipIB,hFL)
    % Second order centered Laplacian  % this term handles viscous
    w = (uFL(ipIB,:,:,:)+uFL(imIB,:,:,:)+uFL(:,ipIB,:,:)+uFL(:,imIB,:,:)+...
    uFL(:,:,ipIB,:)+uFL(:,:,imIB,:)-6*uFL)/(hFL*hFL);
end