function [bif, rads, vpts, froze] = bifurcation_detection(im,mask,tx,ty,maxR)

stepSize = 5;
R = 5:stepSize:maxR;
idx = find(mask==1);
perm_idx = randperm(length(idx));

count = 1;
froze = zeros(size(im));
for i = 1:length(idx)
    ctrp = idx(perm_idx(i));
    if froze(ctrp) == 1
        continue;
    end
    [ctr_r,ctr_c] = ind2sub(size(im),ctrp);
    
    % for debug only
    ctr_r = 192; ctr_c = 425;
%     ctr_r = 102; ctr_c = 282;
    
    for j = 1:length(R)
        radius = R(j);
        [x, y] = circle(ctr_c,ctr_r,radius);
        % get local roi
        r_range = max(1,min(floor(y))-1):min(max(ceil(y))+1,size(im,1));
        c_range = max(1,min(floor(x))-1):min(max(ceil(x))+1,size(im,2));
        roi = im(r_range,c_range);
        [X,Y] = meshgrid(1:size(roi,2),1:size(roi,1));
        profile = interp2(X,Y,roi,x-min(c_range)+1,y-min(r_range)+1);
        % rearrange profile
        [~, pos] = min(profile);
        if pos~=1 && pos~=length(profile) 
            profile = [profile(pos:end); profile(1:pos-1)];
            x = [x(pos:end); x(1:pos-1)];
            y = [y(pos:end); y(1:pos-1)];
        end
        profile = smooth(smooth(profile));
        [xmax,imax,xmin,imin] = extrema(profile);
        imax(xmax<0.1) = [];
        xmax(xmax<0.1) = [];
        if length(imax) < 3
            continue;
        end
        max_pts = [y(imax), x(imax)];
        big_r_range = max(1,r_range(1)-5):min(size(im,1),r_range(end)+5);
        big_c_range = max(1,c_range(1)-5):min(size(im,2),c_range(end)+5);
        bigroi = im(big_r_range,big_c_range);
        roitx = tx(big_r_range,big_c_range);
        roity = ty(big_r_range,big_c_range);
        
        [tx_new1,ty_new1,score1] = ETFstraight(false, 0, roitx, roity,bigroi);
        % *************correct opposite direction**************
        lgc = acos([tx_new1(:),ty_new1(:)]*[0;-1]) < 0.35;
        ty_new1(lgc) = -ty_new1(lgc);
        [tx_new1,ty_new1,score1] = ETFstraight(false, 1,tx_new1,ty_new1,bigroi);
        
        roitx = tx_new1(6:end-5,6:end-5);
        roity = ty_new1(6:end-5,6:end-5);
        
        % change max_pts x,y to roi x,y
        roi_pts = [max_pts(:,1)-min(r_range)+1, max_pts(:,2)-min(c_range)+1];
        roi_pts1 = round(roi_pts);
        pts_direction = zeros(size(roi_pts));
        for num = 1:size(roi_pts,1)
            pts_direction(num,2) = roitx(roi_pts1(num,1),roi_pts1(num,1));
            pts_direction(num,1) = roity(roi_pts1(num,1),roi_pts1(num,1));
        end
        % angle bewteen direction and center-to-point vector
        roi_ctr = [ctr_r-min(r_range)+1, ctr_c-min(c_range)+1];
        vecs = bsxfun(@minus, roi_pts, roi_ctr);
        vecs = bsxfun(@rdivide, vecs, sqrt(sum(vecs.^2,2)));
        angs = acos(sum(pts_direction.*vecs,2))./pi.*180;
        lgc = angs <= 30 | angs >= 150;
        
        % for debug
%         figure
%         imshow(roi,[])
%         hold on
%         quiver(roitx,roity)
%         plot(roi_pts(:,2),roi_pts(:,1),'yo');
        
        if sum(lgc) < 3
            continue;
        else
            tmp1 = max_pts(angs <= 30,:);
            tmp1(:,3) = 1;
            tmp2 = max_pts(angs >= 150,:);
            tmp2(:,3) = -1;
            vpts{count} = [tmp1; tmp2];
            bif(count,:) = [ctr_r,ctr_c];
            rads(count) = R(j);
            [circX, circY] = drawSolidCircle([ctr_r,ctr_c], R(j)+stepSize, size(im));
            circInd = sub2ind(size(im),circY,circX);
            froze(circInd) = 1;
            count = count+1;
            break;
        end
    end
    [circX, circY] = drawSolidCircle([ctr_r,ctr_c], min(R), size(im));
    circInd = sub2ind(size(im),circY,circX);
    froze(circInd) = 1;
    
    if mod(i,100) == 0
        i
        toc
    end
    
end