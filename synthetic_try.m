addpath F:\Dropbox\Code\IPfunctions

% toy experiment
ctr1 = [25,25];
vec1 = [0,-1];
len1 = 15;

ctr2 = [25,25+len1];
vec2 = [-1,-1];
len2 = 12;

ctr3 = [25,25+len1];
vec3 = [-1,1];
len3 = 12;

[x1, y1] = line2D(ctr1, vec1, len1);
[x2, y2] = line2D(ctr2, vec2, len2);
[x3, y3] = line2D(ctr3, vec3, len3);

lgc2 = x2<ctr2(2);
x2(lgc2) = [];
y2(lgc2) = [];

lgc3 = x3<ctr3(2);
x3(lgc3) = [];
y3(lgc3) = [];

% figure
% hold on
% plot(x1,y1,x2,y2,x3,y3)

imgSize = [ceil(max([y1;y2;y3])/10)*10, ceil(max([x1;x2;x3])/10)*10];
mask1 = reconstruct([y1,x1], 4, imgSize);
mask2 = reconstruct([y2,x2], 3, imgSize);
mask3 = reconstruct([y3,x3], 3, imgSize);
mask = mask1 | mask2 | mask3;
% figure
% imshow(mask)
%% add noise
img = double(mask);
img(img==0) = 0.3;
img(img==1) = 0.7;
% img = Norm(imnoise(double(mask),'gaussian',0.5,0.2));
% figure
% imshow(img,[])
%% 
addpath 'helper_function'
addpath F:\Dropbox\Code\experiment

% img = imtophat(1-img,strel('disk',12));
roi = img;
%% example section
% r1 = [344,454,399,484];
% roi1 = img(r1(1):r1(3),r1(2):r1(4));
% roi = roi1;
r2 = [173,399,210,438];
roi2 = img(r2(1):r2(3),r2(2):r2(4));
roi = roi2;
%%
thres = 0.65; sigma = 0; max_L = 15; 
w = 15;
[R_combined, R_theta, rmask] = detect_lines(roi, ones(size(roi)), ones(size(roi)), w, max_L, thres, sigma);

[tx, ty, gmag_norm] = get_tangent(R_combined);
[tx_new, ty_new] = ETF(tx,ty,gmag_norm,roi);
[tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,roi);
[tx_new,ty_new,score1] = ETFstraight(false, 0,tx_new,ty_new,roi);
% *************correct opposite direction**************
lgc = acos([tx_new(:),ty_new(:)]*[0;-1]) < 0.35;
ty_new(lgc) = -ty_new(lgc);
[tx_new,ty_new,score1] = ETFstraight(false, 1,tx_new,ty_new,roi);

[X,Y]=meshgrid(1:size(roi,2),1:size(roi,1));

%%
feat = [roi(:),R_combined(:),tx_new(:),ty_new(:)];%,X(:),Y(:)];
feat(:,1) = Norm(feat(:,1));
feat(:,2) = Norm(feat(:,2));
feat(:,3) = Norm(feat(:,3));
feat(:,4) = Norm(feat(:,4));


[coeff,score] = princomp(feat);
figure
plot(score(:,1),score(:,2),'+');

%%
% mask = zeros(size(roi));
% mask(score(:,1)<0) = 1;
% figure
% imshow(mask);
%%
mask = zeros(size(roi));
mask(score(:,2)>0.2) = 1;
figure
imshow(mask);
%%
mask = zeros(size(roi));
% mask(score(:,2)>0.2 & score(:,1)<0.2 & score(:,1)>-0.4) = 1;
mask(score(:,2)>0.2 & score(:,1)>=0.2) = 1;
% mask(score(:,2)>0.2 & score(:,1)<=-0.4) = 1;
% figure
% imshow(mask);

%% plot vectors
midx = find(mask==1);
[mr,mc] = find(mask==1);
% figure
hold on
quiver(mc,mr,tx_new(midx),ty_new(midx));

%%
feat = [roi(:),R_combined(:),tx_new(:),ty_new(:),X(:),Y(:)];
feat(:,5) = Norm(feat(:,5));
feat(:,6) = Norm(feat(:,6));
[centres, idx] = my_kmeans(feat, 6);
figure
imshow(reshape(idx,size(roi,1),size(roi,2)),[]);