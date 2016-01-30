%% main script
addpath 'helper_function'
 
[I, M, G, rgbI] = openfile('DRIVE', 'training');

for i = 9:9%1:size(I,3)
    % get current instance
    im = I(:,:,i); cmask = M(:,:,i); vmask = G(:,:,i);
    %% preprocessing
    [img, img_raw, pmask] = vessel_preproc(im, cmask);
    img = imtophat(1-img,strel('disk',12));
    
    [tx, ty, gmag_norm] = get_tangent(img);
    
    [tx_new2, ty_new2, tmag1] = ETF_gpu(tx,ty,gmag_norm,img,'etfOriginal');
    [tx_new2, ty_new2, tmag2] = ETF_gpu(tx_new2',ty_new2',gmag_norm,img,'etfOriginal');
    [tx_new2, ty_new2, tmag3] = ETF_gpu(tx_new2',ty_new2',gmag_norm,img,'etfStraight');
    
    tx_new2 = tx_new2';
    ty_new2 = ty_new2';
    
    tmag = sqrt(tx_new2.^2+ty_new2.^2);
    tmag(isnan(tmag)) = 0;
    cc = bwconncomp(tmag);
    bwmask = zeros(size(tmag));
    for j = 1:cc.NumObjects
        idx = cc.PixelIdxList{j};
        if length(idx)>=10
            bwmask(idx) = 1;
        end
    end
    idx = find(bwmask.*cmask);
    tx = tx_new2(idx);
    ty = ty_new2(idx);
    
    [y,x] = ind2sub(size(img),idx);
    figure
    imshow(img.*cmask,[])
    figure
    imshow(img.*cmask,[])
    hold on
    quiver(x,y,tx,ty);
    
%     save(['quiver_results/',num2str(i),'.mat'],'y','x', 'tx','ty');
    
end