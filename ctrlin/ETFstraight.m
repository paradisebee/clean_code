function [r_tx,r_ty,score] = ETFstraight(flag, tf, sm_tx,sm_ty,sm_roi,sm_mask)

% addpath F:\Dropbox\Code\IPfunctions

if nargin < 6
    sm_mask = ones(size(sm_roi));
end
ctrps = find(sm_mask==1);

alpha = 0:pi/4:(pi-pi/4);
[r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, sm_tx, sm_ty, flag);
if tf == 1
[r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, r_tx, r_ty, flag);  
end
% [r_tx,r_ty] = direcCorrect(ctrps,r_tx,r_ty);
end


function [r_tx, r_ty, score] = iter(ctrps, alpha, sm_roi, sm_tx, sm_ty, cross)

[r,c] = ind2sub(size(sm_roi),ctrps);
r_tx = zeros(size(sm_tx));
r_ty = zeros(size(sm_ty));
score = r_tx;

len = 4;

[tx, ty, gmag_norm] = get_tangent(sm_roi);

for j = 1:length(ctrps)
    
    ctrp = [r(j), c(j)];
    vec = [-sm_ty(ctrp(1),ctrp(2)),sm_tx(ctrp(1),ctrp(2))];
    [x, y] = line2D(ctrp, vec, len);
    x = floor(x+0.5); y = floor(y+0.5);
    lgc = y<=0 | x<=0 | y>size(sm_roi,1) | x>size(sm_roi,2);
    y(lgc) = [];
    x(lgc) = [];
    
    idx = sub2ind(size(sm_roi),y,x);
    idx = unique(idx);
    if length(idx)==1
        continue;
    end

    cand = zeros(size(sm_roi));
    cand(idx) = 1;
    
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
    
    
    if cross==true
    %% cross
    vec2 = [vec(2), -vec(1)];
    [x2, y2] = line2D(ctrp, vec2, len);
    x2 = floor(x2+0.5); y2 = floor(y2+0.5);
    lgc = y2<=0 | x2<=0 | y2>size(sm_roi,1) | x2>size(sm_roi,2);
    y2(lgc) = [];
    x2(lgc) = [];
    

    sub = [y2,x2];
    idx = sub2ind(size(sm_roi),sub(:,1),sub(:,2));
    idx_diff = idx(2:end)-idx(1:end-1);
    loc = find(idx_diff==0)+1;
    
    sub(loc,:) = [];
    idx(loc) = [];
    ctr_pos = find(sub(:,1)==ctrp(1) & sub(:,2)==ctrp(2));
    
    
    
    x3 = x2(1:ctr_pos);
    y3 = y2(1:ctr_pos); 
    x4 = x2(ctr_pos+1:end);
    y4 = y2(ctr_pos+1:end);
    
    % conditions
    inten = sm_roi(idx);
    lg1 = inten<0.2;
    % lg1_list = bwconncomp1D(lg1);
    lg2 = inten<0.4;
    % lg2_list = bwconncomp1D(lg2);
    lg3 = inten<0.5;
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
                fx = -vec2(1); fy = -vec(2);
            else
                fx = vec2(1); fy = vec2(2);
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
                fx = -vec2(1); fy = -vec2(2);
            else
                fx = vec2(1); fy = vec2(2);
            end
            fx = fx*w1; fy = fy*w1; s = 2;
        elseif sum(lg1)
            lg1pos = find(lg1);
            [~,ii] = min(abs(ctr_pos-lg1pos));
            if sum(y4==sub(lg1pos(ii),1) & x4==sub(lg1pos(ii),2))
                fx = -vec2(1); fy = -vec2(2);
            else
                fx = vec2(1); fy = vec2(2);
            end
            fx = fx*w2; fy = fy*w2; s = 3;
        end
        flag = 1;
    else
        if sum(lg1) == 0 && sum(lg2) == 0 && sum(lg3) == 0
            flag = 0;
            if min(sm_roi(cand==1))<0.5
                fx = 0; fy = 0; s = 1;
                flag = 1;
            end
        elseif sum(lg3)
            lg3pos = find(lg3);
            [~,ii] = min(abs(ctr_pos-lg3pos));
            if sum(y4==sub(lg3pos(ii),1) & x4==sub(lg3pos(ii),2))
                fx = -vec2(1); fy = -vec2(2);
            else
                fx = vec2(1); fy = vec2(2);
            end
            fx = fx*w1; fy = fy*w1; flag = 1; s = 2;
        elseif sum(lg2)
            lg2pos = find(lg2);
            [~,ii] = min(abs(ctr_pos-lg2pos));
            if sum(y4==sub(lg2pos(ii),1) & x4==sub(lg2pos(ii),2))
                fx = -vec2(1); fy = -vec2(2);
            else
                fx = vec2(1); fy = vec2(2);
            end
            fx = fx*w2; fy = fy*w2; flag = 1; s = 3;
        elseif sum(lg1)
            lg1pos = find(lg1);
            [~,ii] = min(abs(ctr_pos-lg1pos));
            if sum(y4==sub(lg1pos(ii),1) & x4==sub(lg1pos(ii),2))
                fx = -vec2(1); fy = -vec2(2);
            else
                fx = vec2(1); fy = vec2(2);
            end
            fx = fx; fy = fy; flag = 1; s = 4;
        end
    end
    
    if flag == 1
        r_tx(ctrp(1),ctrp(2)) = r_tx(ctrp(1),ctrp(2))+fx;
        r_ty(ctrp(1),ctrp(2)) = r_ty(ctrp(1),ctrp(2))+fy;
    else
        r_tx(ctrp(1),ctrp(2)) = 0;
        r_ty(ctrp(1),ctrp(2)) = 0;
    end
    score(ctrp(1),ctrp(2)) = s;
    end
    
end

r_mag = sqrt(r_tx.^2+r_ty.^2);
r_tx = r_tx./sqrt(r_tx.^2+r_ty.^2);
r_ty = r_ty./sqrt(r_tx.^2+r_ty.^2);
r_tx(r_mag==0) = 0;
r_ty(r_mag==0) = 0;

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

%%
function [tx,ty] = direcCorrect(ctrps,tx,ty)

[r,c] = ind2sub(size(tx),ctrps);
for j = 1:length(ctrps)
    
    ctrp = [r(j), c(j)];
    vec = [tx(ctrp(1),ctrp(2)),ty(ctrp(1),ctrp(2))];
    [x, y] = line2D(ctrp, vec, 4);
    x = floor(x+0.5); y = floor(y+0.5);
    lgc = y<=0 | x<=0 | y>size(tx,1) | x>size(tx,2);
    y(lgc) = [];
    x(lgc) = [];
    
    idx = sub2ind(size(tx),y,x);
    idx = unique(idx);
    if length(idx)==1
        continue;
    end

    win_tx = tx(idx);
    win_ty = ty(idx);
    ang = vec*[win_tx,win_ty]';
    if sum(ang<0)>4
        tx(ctrp(1),ctrp(2)) = -tx(ctrp(1),ctrp(2));
        ty(ctrp(1),ctrp(2)) = -ty(ctrp(1),ctrp(2));
    end
end

end