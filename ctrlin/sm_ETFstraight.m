function [ctr_tx, ctr_ty] = sm_ETFstraight(im, tx, ty, ctrp, win_sz)

r_range = ctrp(1)-win_sz:ctrp(1)+win_sz;
c_range = ctrp(2)-win_sz:ctrp(2)+win_sz;
roi = im(r_range,c_range);
roi_tx = tx(r_range,c_range);
roi_ty = ty(r_range,c_range);
roi_ctrp = [ctrp(1)-min(r_range)+1, ctrp(2)-min(c_range)+1];

[tx_new1,ty_new1,score1] = ETFstraight(false, 0, roi_tx, roi_ty,roi);

% *************correct opposite direction**************
lgc = acos([tx_new1(:),ty_new1(:)]*[0;-1]) < 0.35;
ty_new1(lgc) = -ty_new1(lgc);
[tx_new1,ty_new1,score1] = ETFstraight(false, 1,tx_new1,ty_new1,roi);

ctr_tx = tx_new1(roi_ctrp(1),roi_ctrp(2));
ctr_ty = ty_new1(roi_ctrp(1),roi_ctrp(2));