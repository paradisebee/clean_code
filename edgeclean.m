function [outIm, edges] = edgeclean(edgeIm, len)

bw = double(bwmorph(edgeIm, 'thin'));
num_of_edge_points = sum(bw(:));
search_map = bw;
pidx = find(search_map==1);
cp = pidx(round(1+(num_of_edge_points-1).*rand(1)));

% [tx, ty, gmag_norm] = get_tangent(bw);

% set kernel and search map
% ----------- one-side edge ---------------
knl1 = [1,2,3;5,8,13;21,34,55];
knl_map1 = bsxfun(@plus, knl1(:), knl1(:)');
knl_map1(:,[1:4,6:9]) = 0;
% ----------- two-side edge ----------------
knl2 = [1,2,3;5,0,13;21,34,55];
knl_map2 = bsxfun(@plus, knl2(:), knl2(:)');
knl_map2 = knl_map2+8;
for i = 1:9
    knl_map2(i,1:i) = 0;
end
knl_map2(:,5) = 0; knl_map2(5,:) = 0;
max_knl2 = max(knl_map2(:));


% get tangent vector of edge points
il = imfilter(bw, knl1, 'same');
il(bw~=1) = 0;
idx = find(il~=0);
[ir,ic] = find(il~=0);

tx = zeros(size(bw)); ty = zeros(size(bw));
neigh = [-1,-1;
    0,-1;
    1,-1;
    -1, 0;
    0, 0;
    1, 0;
    -1, 1;
    0, 1;
    1, 1;];
for i = 1:length(idx)
    val = il(idx(i));
    [r,c] = find(knl_map2==val);
    if isempty(r)
        [r,c] = find(knl_map1==val);
        if ~isempty(r)
            pxy = [neigh(r,2), neigh(r,1)];
        else
            continue;
        end
    end
    v1 = [neigh(r,2),neigh(r,1)];
    v2 = [neigh(c,2),neigh(c,1)];
    ang12 = acos(v1*v2'/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2))));
    if ang12 > pi/2
        v1 = -v1;
    end
    pxy = v1 + v2;
    tx(idx(i)) = pxy(1); 
    ty(idx(i)) = pxy(2);
end
tx = tx./sqrt(tx.^2+ty.^2+1e-5);
ty = ty./sqrt(tx.^2+ty.^2+1e-5);
% 
% figure
% imshow(bw)
% hold on
% quiver(ic,ir,tx(idx),ty(idx));

%% main cleaning
rnd_idx = randperm(num_of_edge_points);
npidx = pidx(rnd_idx);
froze_map = zeros(size(bw));
EDGE = {};
for i = 1:length(rnd_idx)
    curP = npidx(i);
    if froze_map(curP) == 1
        continue;
    end
    froze_map(curP) = 1;
    curT = [tx(curP), ty(curP)];
    sub_edge = curP;
    while length(sub_edge) < len
        [cr,cc] = ind2sub(size(bw),curP);
        next_sub = round([cr+curT(2), cc+curT(1)]);
        nextP = sub2ind(size(bw),next_sub(1),next_sub(2));
        if bw(nextP)==1 && ~sum(sub_edge==nextP) && froze_map(nextP)~=1
            sub_edge = [sub_edge,nextP];
            curP = nextP;
            curT = [tx(curP), ty(curP)];
        else
            EDGE = [EDGE,sub_edge];
            break;
        end
    end
    
    
end

figure
imshow(edgeIm,[])
hold on
for i = 1:length(EDGE)
    s = EDGE{i};
    [sr,sc] = ind2sub(size(bw),s);
    plot(sc,sr,'xr');
%     pause(0.5);
end

figure;
