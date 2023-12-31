function [my_vertices,my_faces] = makeWallPlotVars(wallR,Aeff)

minwx = min(wallR(:,1))-Aeff;
maxwx = max(wallR(:,1))+Aeff;
minwy = min(wallR(:,2))-Aeff;
maxwy = max(wallR(:,2))+Aeff;
minwz = min(wallR(:,3))-Aeff;
maxwz = max(wallR(:,3))+Aeff;
my_vertices = [minwx minwy minwz;...
               minwx maxwy minwz;...
               maxwx maxwy minwz;...
               maxwx minwy minwz;...
               minwx minwy maxwz;...
               minwx maxwy maxwz;...
               maxwx maxwy maxwz;...
               maxwx minwy maxwz];
my_faces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
end