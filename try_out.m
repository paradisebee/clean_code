addpath C:\Users\Bichao\Desktop\clean_code\helper_function

load edge.mat
af = Norm(af);
Ieb = Norm(Ieb);
figure
imshow(af,[])

% use heaviside function 
epsilon = 0.3;
h_thres = 0.8;
heav_f = @(a,h_thres) 0.5*(1+2/pi*atan((a-h_thres)/epsilon));

hf = heav_f(Ieb,h_thres);

% fh = @(a,b) (a==b)*0.2+(a~=b)*1;
h_thres = 0.2;
fh = @(a,b) ((heav_f(abs(a-b),h_thres)>0.5)*max(heav_f(abs(a-b),h_thres),1)+...
    (1-heav_f(abs(a-b),h_thres)>0.5)*min(heav_f(abs(a-b),h_thres),0.01));
pairwise = setSmoothness2Dfunc(hf,0.01,0.5,'8nbh',fh);

% unary = zeros(2,length(af(:)));
% unary(1,:) = 10*(af(:))';
% unary(2,:) = 10*(1-af(:))';

prob_obj = zeros(size(af));
% prob_obj(af>0.4) = 0.9;
% prob_obj(af<=0.4) = 0.1;
prob_obj = af;

prob_bkg = zeros(size(af));
% prob_bkg(af>0.4) = 0;
% prob_bkg(af<=0.4) = 1;
prob_bkg = 1-af;

unary = zeros(2,length(af(:)));
unary(1,:) = -log(prob_obj(:));
unary(2,:) = -log(prob_bkg(:));
unary(unary==0) = 1e-5;
unary(unary==inf) = 100;
% unary(unary~=0.5) = 0.5;


inp = zeros(length(af(:)),1);
inp = 1-double(af>0.4);
labelcost = ones(2) - eye(2);

[labels E Eafter] = GCMex(inp(:), single(unary), pairwise, single(labelcost),0);%segclass
gmask = reshape(labels, size(af,1), size(af,2));
figure; 
imshow(1-gmask);