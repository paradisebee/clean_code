function [r_tx,r_ty,score] = ETFCross(sm_tx,sm_ty,sm_roi,sm_mask)
% close all
% tic
% load im.mat
% sm_roi = roi_corr(346:458, 401:478);
% sm_mask = roi_mask(346:458, 401:478);
% sm_tx = tx_new(346:458, 401:478);
% sm_ty = ty_new(346:458, 401:478);

% sm_roi = roi_corr(439:508, 325:380);
% sm_mask = roi_mask(439:508, 325:380);
% sm_tx = tx_new(439:508, 325:380);
% sm_ty = ty_new(439:508, 325:380);


% figure
% imshow(sm_roi,[]);


% alpha = 0:pi/12:(pi-pi/12);
alpha = 0:pi/4:(pi-pi/4);
% alpha = alpha./pi.*180;



if nargin < 4
    sm_mask = ones(size(sm_roi));
end
ctrps = find(sm_mask==1);



[r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, sm_tx, sm_ty);
[r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, r_tx, r_ty);    

end

function [r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, sm_tx, sm_ty)

[r,c] = ind2sub(size(sm_roi),ctrps);
r_tx = zeros(size(sm_tx));
r_ty = zeros(size(sm_ty));
score = r_tx;

[tx, ty, gmag_norm] = get_tangent(sm_roi);

for j = 1:length(ctrps)

% ctrp = [41,43];
% ctrp = [11,40];
% ctrp = [5, 56];
% ctrp = [59,21];
% ctrp = [78,60];
% ctrp = [37,70];
% ctrp = [8,55];
ctrp = [r(j), c(j)];

similarity = [];

for i = 1:length(alpha)
theta = alpha(i);
% ctrp = round(size(sm_roi)/2);
radius = 0:4;

y1 = round(ctrp(1)+radius*sin(theta));
x1 = round(ctrp(2)+radius*cos(theta));
lgc1 = y1<=0 | x1<=0 | y1>size(sm_roi,1) | x1>size(sm_roi,2);
y1(lgc1) = [];
x1(lgc1) = [];

y2 = round(ctrp(1)-radius*sin(theta));
x2 = round(ctrp(2)-radius*cos(theta));
lgc2 = y2<=0 | x2<=0 | y2>size(sm_roi,1) | x2>size(sm_roi,2);
y2(lgc2) = [];
x2(lgc2) = [];


idx = sub2ind(size(sm_roi),[y1,y2],[x1,x2]);
idx = unique(idx);


se = strel('disk',1);

cand = zeros(size(sm_roi));
cand(idx) = 1;
% cand2 = imdilate(imdilate(imdilate(cand,se),se),se);
% cand2 = imdilate(imdilate(cand,se),se);
cand2 = imdilate(cand,se);
% cand2 = cand;

% figure
% imshow(cand);
% figure
% imshow(cand2);

cx = sm_tx(cand2==1);
cy = sm_ty(cand2==1);

angles = acos((cx*cos(theta)+cy*sin(theta))./(sqrt(cx.^2+cy.^2)))./pi.*180;
l1 = angles;
l1(angles>135) = l1(angles>135)-180;
l1 = l1./180.*pi;

l2 = angles;
if (sum((l2>=45) & (l2<=135)) >= length(angles)/2) || length(idx) == 1
    similarity(i) = 1;
else
    similarity(i) = var(l1)*(1-length(idx)/10);
end


end

[minsim, pos] = min(similarity);
theta = alpha(pos);

y1 = round(ctrp(1)+radius*sin(theta));
x1 = round(ctrp(2)+radius*cos(theta));
lgc1 = y1<=0 | x1<=0 | y1>size(sm_roi,1) | x1>size(sm_roi,2);
y1(lgc1) = [];
x1(lgc1) = [];

y2 = round(ctrp(1)-radius*sin(theta));
x2 = round(ctrp(2)-radius*cos(theta));
lgc2 = y2<=0 | x2<=0 | y2>size(sm_roi,1) | x2>size(sm_roi,2);
y2(lgc2) = [];
x2(lgc2) = [];


idx = sub2ind(size(sm_roi),[y1,y2],[x1,x2]);
idx = unique(idx);
if length(idx)==1
    continue;
end

cand = zeros(size(sm_roi));
cand(idx) = 1;
% cand2 = imdilate(imdilate(cand,se),se);
cand2 = imdilate(cand,se);
% cand2 = cand;
% idx = find(cand2==1);


%%

win_gmag = gmag_norm(idx);
win_tx = sm_tx(idx);
win_ty = sm_ty(idx);
t_cur = [win_tx(:),win_ty(:)];

wm = mag_weight(gmag_norm(ctrp(1),ctrp(2)),win_gmag(:),1);
wd = dire_weight([sm_tx(ctrp(1),ctrp(2)),sm_ty(ctrp(1),ctrp(2))],[win_tx(:),win_ty(:)]);

temp = sum(bsxfun(@times, t_cur, wm.*wd));
if sqrt(sum(temp.^2))~=0
    temp = temp./sqrt(sum(temp.^2));
end

if sum(wd)/9 < -0.1
    temp = -temp;
end

r_tx(ctrp(1),ctrp(2)) = temp(1);
r_ty(ctrp(1),ctrp(2)) = temp(2);


% figure
% imshow(cand2);
% hold on
% quiver(ctrp(2),ctrp(1),temp(1),temp(2));

%% cross
y3 = round(ctrp(1)+radius*sin((theta/pi*180+90)/180*pi));
x3 = round(ctrp(2)+radius*cos((theta/pi*180+90)/180*pi));
lgc3 = y3<=0 | x3<=0 | y3>size(sm_roi,1) | x3>size(sm_roi,2);
y3(lgc3) = [];
x3(lgc3) = [];

y4 = fliplr(round(ctrp(1)-radius*sin((theta/pi*180+90)/180*pi)));
x4 = fliplr(round(ctrp(2)-radius*cos((theta/pi*180+90)/180*pi)));
lgc4 = y4<=0 | x4<=0 | y4>size(sm_roi,1) | x4>size(sm_roi,2);
y4(lgc4) = [];
x4(lgc4) = [];

sub = [[y4,y3]',[x4,x3]'];
idx = sub2ind(size(sm_roi),sub(:,1),sub(:,2));
idx_diff = idx(2:end)-idx(1:end-1);
loc = find(idx_diff==0)+1;

sub(loc,:) = [];
idx(loc) = [];
ctr_pos = find(sub(:,1)==ctrp(1) & sub(:,2)==ctrp(2));

% conditions
inten = sm_roi(idx);
lg1 = inten>=0.4;
% lg1_list = bwconncomp1D(lg1);
lg2 = inten>=0.25;
% lg2_list = bwconncomp1D(lg2);
lg3 = inten>=0.15;
% lg3_list = bwconncomp1D(lg3);

w1 = 0.3; w2 = 0.6; s = 0;
% w1 = -0.8; w2 = -1;
if lg1(ctr_pos) == 1
    fx = 0; fy = 0; s = 1;
    flag = 1;
elseif lg2(ctr_pos) == 1
    if sum(lg1) == 0
        fx = 0; fy = 0; s = 1;
    else
        lg1pos = find(lg1);
        [~,ii] = min(abs(ctr_pos-lg1pos));
        if sum(y4==sub(lg1pos(ii),1) & x4==sub(lg1pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx*w1; fy = fy*w1; s = 2;
    end
    flag = 1;
elseif lg3(ctr_pos) == 1
    if sum(lg1) == 0 && sum(lg2) == 0
        fx = 0; fy = 0; s = 1;
    elseif sum(lg2)
        lg2pos = find(lg2);
        [~,ii] = min(abs(ctr_pos-lg2pos));
        if sum(y4==sub(lg2pos(ii),1) & x4==sub(lg2pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx*w1; fy = fy*w1; s = 2;
    elseif sum(lg1)
        lg1pos = find(lg1);
        [~,ii] = min(abs(ctr_pos-lg1pos));
        if sum(y4==sub(lg1pos(ii),1) & x4==sub(lg1pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx*w2; fy = fy*w2; s = 3;
    end
    flag = 1;
else
    if sum(lg1) == 0 && sum(lg2) == 0 && sum(lg3) == 0
        flag = 0;
        if minsim<=0.1 || min(sm_roi(cand==1))<0.5
            fx = 0; fy = 0; s = 1;
            flag = 1;
        end
    elseif sum(lg3)
        lg3pos = find(lg3);
        [~,ii] = min(abs(ctr_pos-lg3pos));
        if sum(y4==sub(lg3pos(ii),1) & x4==sub(lg3pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx*w1; fy = fy*w1; flag = 1; s = 2;
    elseif sum(lg2)
        lg2pos = find(lg2);
        [~,ii] = min(abs(ctr_pos-lg2pos));
        if sum(y4==sub(lg2pos(ii),1) & x4==sub(lg2pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx*w2; fy = fy*w2; flag = 1; s = 3;
    elseif sum(lg1)
        lg1pos = find(lg1);
        [~,ii] = min(abs(ctr_pos-lg1pos));
        if sum(y4==sub(lg1pos(ii),1) & x4==sub(lg1pos(ii),2))
            fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
        else
            fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
        end
        fx = fx; fy = fy; flag = 1; s = 4;
    end
end
% flag = 0;
% if lg1(ctr_pos)==1
%     fx = 0; fy = 0;
%     flag = 1;
% elseif lg2(ctr_pos)==1
%     lg1pos = find(lg1);
%     if isempty(lg1pos)
%         fx = 0; fy = 0;
%     else
%         [~,ii] = min(abs(ctr_pos-lg1pos));
%         if sum(y4==sub(ii,1) & x4==sub(ii,2))
%             fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
%         else
%             fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
%         end
%         fx = fx*0.3; fy = fy*0.3;
%     end
%     flag = 1;
% elseif lg3(ctr_pos)==1
%     lg1pos = find(lg1);
%     lg2pos = find(lg2);
%     if isempty(lg1pos) || isempty(lg2pos)
%         flag = 0;
%     elseif sum(lg1pos)
%         [~,ii] = min(abs(ctr_pos-lg1pos));
%         if sum(y4==sub(ii,1) & x4==sub(ii,2))
%             fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
%         else
%             fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
%         end
%         fx = fx; fy = fy;
%         flag = 1;
%     else
%         [~,ii] = min(abs(ctr_pos-lg2pos));
%         if sum(y4==sub(ii,1) & x4==sub(ii,2))
%             fx = cos((theta/pi*180+90)/180*pi); fy = sin((theta/pi*180+90)/180*pi);
%         else
%             fx = -cos((theta/pi*180+90)/180*pi); fy = -sin((theta/pi*180+90)/180*pi);
%         end
%         fx = fx*0.6; fy = fy*0.6;
%         flag = 1;
%     end
% else 
%     if 
% end

if flag == 1
    r_tx(ctrp(1),ctrp(2)) = r_tx(ctrp(1),ctrp(2))+fx;
    r_ty(ctrp(1),ctrp(2)) = r_ty(ctrp(1),ctrp(2))+fy;
else
    r_tx(ctrp(1),ctrp(2)) = 0;
    r_ty(ctrp(1),ctrp(2)) = 0;
end
score(ctrp(1),ctrp(2)) = s;

end

r_mag = sqrt(r_tx.^2+r_ty.^2);
r_tx = r_tx./sqrt(r_tx.^2+r_ty.^2);
r_ty = r_ty./sqrt(r_tx.^2+r_ty.^2);
r_tx(r_mag==0) = 0;
r_ty(r_mag==0) = 0;

% figure
% imshow(sm_roi,[]);
% hold on
% quiver(sm_tx,sm_ty);
% figure
% imshow(sm_roi,[]);
% hold on
% quiver(r_tx,r_ty);

% toc
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

