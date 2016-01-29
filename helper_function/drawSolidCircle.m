function [circX, circY] = drawSolidCircle(ctr, radius, imgSize)

[X,Y] = meshgrid(ctr(2)-radius:ctr(2)+radius, ctr(1)-radius:ctr(1)+radius);
lgc = (X-ctr(2)).^2+(Y-ctr(1)).^2<radius^2 & X>=1 & X<=imgSize(2) & Y>=1 & Y<=imgSize(1);
circX = X(lgc);
circY = Y(lgc);