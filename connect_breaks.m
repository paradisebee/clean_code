addpath 'helper_function'

[I, M, G, rgbI] = openfile('DRIVE', 'test');

%% set parameters for line detection 
thres = 0.65; sigma = 0; max_L = 15; 
w = 15;

h = figure;
for i = 1:size(I,3)
    % get current instance
    im = I(:,:,i); cmask = M(:,:,i); vmask = G(:,:,i);
    %% preprocessing
    [img, img_raw, pmask] = vessel_preproc(im, cmask);
    img = imtophat(1-img,strel('disk',12)); 
    % line detection
    [R_combined, R_theta, rmask] = detect_lines(img, pmask, cmask, w, max_L, thres, sigma);
    
    
    [tx, ty, gmag_norm] = get_tangent(R_combined);
    [tx_new, ty_new] = ETF(tx,ty,gmag_norm,img);
    [tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,img);
%     k = gca(h);
%     imshow(R_combined,[])
%     alphamask(vmask, [0,1,0], 0.1, k);
%     
%     pause;
    
    
    %% exclud pixels which line detection response are less than 0
    R1 = R_combined;
    R1(R_combined<0) = 0;
    
    %% exclude pixels which direction vectors are not consistant
    [r,c] = size(im);
    [X, Y] = meshgrid(1:c,1:r);
    idx = find(R1.*pmask~=0);
    ep_mask = zeros(size(im));
    tic
    for j = 1:length(idx)
        ctx = tx_new(idx(j)); cty = ty_new(idx(j));
        [cr,cc] = ind2sub(size(R1),idx(j));
        [x, y, pos] = line2D([cr,cc], [ctx,cty], 4);
        htx=interp2(X,Y,tx,x,y);
        hty=interp2(X,Y,ty,x,y);
        htx = htx./sqrt(htx.^2+hty.^2);
        hty = hty./sqrt(htx.^2+hty.^2);
        ang_dif = acos([ctx,cty]*[htx,hty]')./pi.*180;
        ang_dif(isnan(ang_dif)) = 150;
        if sum(ang_dif<45) > 5
            ep_mask(idx(j)) = 1;
        end
        if mod(j,1000)==0
            toc
        end
    end
    
    pause;
end