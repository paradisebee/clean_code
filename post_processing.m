function output = post_processing(img, result, pmask)

% the threshold is the minimum value of 12% highest values in img
inten = sort(img(pmask==1),'descend');
thres = inten(round(length(inten)*0.10));
%
im_obj = get_small_vessels(result, pmask);
se = strel('disk',1);
im_obj2 = imerode(im_obj,se);
bdy = im_obj-im_obj2;
pos_idx = img(find(bdy==1));
IDX = find(bdy);
output = result; 
output(IDX(pos_idx<thres))=0;