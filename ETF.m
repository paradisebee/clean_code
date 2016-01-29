%% ETF
function [tx_new, ty_new, tmagnew] = ETF(tx,ty,gmag_norm,im)
win_sz = 3;

tx_new = zeros(size(tx));
ty_new = zeros(size(ty));

[c,r] = meshgrid(1:size(im,2),1:size(im,1));
for i = 1:length(c(:))
    cmin = c(i)-floor(win_sz/2);
    if cmin <= 0
        cmin = 1;
    end
    cmax = c(i)+floor(win_sz/2);
    if cmax > size(im,2)
        cmax = size(im,2);
    end
    rmin = r(i)-floor(win_sz/2);
    if rmin <= 0
        rmin = 1;
    end
    rmax = r(i)+floor(win_sz/2);
    if rmax > size(im,1)
        rmax = size(im,1);
    end
    
    win_gmag = gmag_norm(rmin:rmax,cmin:cmax);
    win_tx = tx(rmin:rmax,cmin:cmax);
    win_ty = ty(rmin:rmax,cmin:cmax);
    t_cur = [win_tx(:),win_ty(:)];
    
    wm = mag_weight(gmag_norm(r(i),c(i)),win_gmag(:),1);
    wd = dire_weight([tx(r(i),c(i)),ty(r(i),c(i))],[win_tx(:),win_ty(:)]);
    
    temp = sum(bsxfun(@times, t_cur, wm.*wd));
    tx_new(r(i),c(i)) = temp(1);
    ty_new(r(i),c(i)) = temp(2);
%     tx_new(r(i),c(i)) = sum(wm);
%     ty_new(r(i),c(i)) = sum(wd);
end

% normalize new tx ty
tmagnew = sqrt(tx_new.^2+ty_new.^2);
tx_new = tx_new./tmagnew;
ty_new = ty_new./tmagnew;
tx_new(tmagnew==0) = 0;
ty_new(tmagnew==0) = 0;

end

%%
function wm = mag_weight(x,y,eta)

% wm = 0.5.*(1+tanh(eta.*(y-x)));
wm = 0.5.*(1+y-x);

end

%% 
function wd = dire_weight(x,y)
% combine phi and wd in paper
wd = y*(x');

end