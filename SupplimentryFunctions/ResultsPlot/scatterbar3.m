% EXAMPLE 1:
close all
clc
clear

y=[1 2 3 1 2 3 1 2 3];
x=[1 1 1 2 2 2 3 3 3];
z=[-1 2 3 6 -5 4 7 8 9];

MaxColor = max(z);

scatterbar3_test(x, y, z, 1, 1, MaxColor)


function scatterbar3_test(X,Y,Z,width_x, width_y, color)
%SCATTERBAR3   3-D scatter bar graph.
%   SCATTERBAR3(X,Y,Z,WIDTH) draws 3-D bars of height Z at locations X and Y with width WIDTH.
%
%   X, Y and Z must be of equal size.  If they are vectors, than bars are placed
%   in the same fashion as the SCATTER3 or PLOT3 functions.
%
%   If they are matrices, then bars are placed in the same fashion as the SURF
%   and MESH functions.
%
%   The colors of each bar read from the figure's colormap according to the bar's height.
%
%   NOTE:  For best results, you should use the 'zbuffer' renderer.  To set the current
%   figure renderer to 'zbuffer' use the following command:
%
   set(gcf,'renderer','zbuffer')

%
%    % EXAMPLE 2:
%    [X,Y]=meshgrid(-1:0.25:1);
%    Z=2-(X.^2+Y.^2);
%    scatterbar3(X,Y,Z,0.2)
%    colormap(hsv)
%
%    % EXAMPLE 3:
%    t=0:0.1:(2*pi);
%    x=cos(t);
%    y=sin(t);
%    z=sin(t);
%    scatterbar3(x,y,z,0.07)

% By Mickey Stahl - 2/25/02
% Engineering Development Group
% Aspiring Developer

[r,c]=size(Z);
for j=1:r
    for k=1:c
        if ~isnan(Z(j,k))
            drawbar(X(j,k),Y(j,k),Z(j,k),width_x/2, width_y/2, color)
        end
    end
end

zlim=[min(Z(:)) max(Z(:))];
if zlim(1)>0,zlim(1)=0;end
if zlim(2)<0,zlim(2)=0;end
axis([min(X(:))-width_x max(X(:))+width_x min(Y(:))-width_y max(Y(:))+width_y zlim])
caxis([min(Z(:)) max(Z(:))])
end

function drawbar(x,y,z,width_x, width_y, MaxColor)

h(1)=patch([-width_x -width_x width_x width_x]+x,[-width_y width_y width_y -width_y]+y, [0 0 0 0], 'b', 'EdgeAlpha', 0.5);
h(2)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[-1 -1 -1 -1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', 0.5);
h(3)=patch(width_x.*[-1 -1 -1 -1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', 0.5);
h(4)=patch([-width_x -width_x width_x width_x]+x, [-width_y width_y width_y -width_y]+y, [z z z z], 'b', 'EdgeAlpha', 0.5);
h(5)=patch(width_x.*[-1 -1 1 1]+x, width_y.*[1 1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', 0.5);
h(6)=patch(width_x.*[1 1 1 1]+x, width_y.*[-1 -1 1 1]+y, z.*[0 1 1 0], 'b', 'EdgeAlpha', 0.5);

if z > 0
    block_color = [1, 1- z/MaxColor,...
                      1- z/MaxColor];
else
    block_color = [1 + z/MaxColor,...
                   1 + z/MaxColor, 1];
end
set(h, 'facecolor', block_color)
end