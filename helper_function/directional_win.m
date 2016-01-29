function [him,htx,hty,dim] = directional_win(im, tx, ty, ctr, width, height)

% width: half length of the window size perpendicular to the direction
%        vector
% height: half length of the window size along the directional vector

ctr_tx = tx(ctr(1),ctr(2));
ctr_ty = ty(ctr(1),ctr(2));

[mwx, mwy, pos] = line2D(ctr, [ctr_tx, ctr_ty], width);

hx = []; hy = [];
for i = 1:length(mwx)
    [mhx, mhy] = line2D([mwy(i),mwx(i)], [-ctr_ty, ctr_tx], height);
    hx = [hx, mhx];
    hy = [hy, mhy];
end

[r,c] = size(im);
[X, Y] = meshgrid(1:c,1:r);

him=interp2(X,Y,im,hx,hy);
htx=interp2(X,Y,tx,hx,hy);
hty=interp2(X,Y,ty,hx,hy);
ang_dif = [0,-1]-[ctr_tx,ctr_ty];
htx = htx+ang_dif(1);
hty = hty+ang_dif(2);
htx = htx./sqrt(htx.^2+hty.^2);
hty = hty./sqrt(htx.^2+hty.^2);
loc = isnan(him);
him(loc) = 0; 
htx(loc) = 0;
hty(loc) = 0;
figure
imshow(him,[])
hold on
quiver(htx,hty)

dim = [hx(1,1),hy(1,1);hx(1,end),hy(1,end);
       hx(end,1),hy(end,1);hx(end,end),hy(end,end)];
   
figure
imshow(im,[])
hold on
hold on
plot([dim(1,1),dim(3,1)],[dim(1,2),dim(3,2)],'r')
plot([dim(1,1),dim(2,1)],[dim(1,2),dim(2,2)],'r')
plot([dim(3,1),dim(4,1)],[dim(3,2),dim(4,2)],'r')
plot([dim(2,1),dim(4,1)],[dim(2,2),dim(4,2)],'r')