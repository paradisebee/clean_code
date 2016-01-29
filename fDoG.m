function Hmap = fDoG(sigma_C, sigma_M, im, tx, ty, mask)

% addpath F:\Dropbox\Code\IPfunctions
sigma_S = 1.6*sigma_C;
rho = 0.99;
T = 3*sigma_C;
t = (-T:T)';
Gc = 1/(sqrt(2*pi)*sigma_C)*exp(-t.^2/(2*sigma_C^2));
Gs = 1/(sqrt(2*pi)*sigma_S)*exp(-t.^2/(2*sigma_S^2));
ft = Gc-rho.*Gs;
S = 3*sigma_M;
s = (-S:S)';
Spos = find(s==0);
Gm = 1/(sqrt(2*pi)*sigma_M)*exp(-s.^2/(2*sigma_M^2));
[r,c] = find(mask==1);

Hmap = zeros(size(im));
for i = 1:length(r)
    ctr = [r(i), c(i)];
    vec = [tx(ctr(1),ctr(2)), ty(ctr(1),ctr(2))];
    if vec(1)==0 && vec(2)==0
        continue;
    end
    [lx, ly, ls_pos] = get_parallel_line(ctr, tx, ty, S);
    h = zeros(length(lx),1);
    for j = 1:length(lx)
        cur_pt = [ly(j),lx(j)];
        vec = [tx(cur_pt(1),cur_pt(2)), ty(cur_pt(1),cur_pt(2))];
        [x, y, pos] = line2D(cur_pt, vec, T);
        val = zeros(size(ft));
        lgc = round(y) >= 1 & round(y)<=size(im,1) & ...
                round(x) >= 1 & round(x)<=size(im,2);
        h(j) = sum(ft(lgc).*im(sub2ind(size(im), round(y(lgc)),round(x(lgc)))));
    end
    if length(h) ~= length(Gm)
        Hmap(ctr(1),ctr(2)) = sum(Gm(Spos-ls_pos+1:Spos+length(h)-ls_pos).*h);
    else
        Hmap(ctr(1),ctr(2)) = sum(Gm.*h);
    end
end


% miscellaneous
function [x, y, pos] = get_parallel_line(ctr, tx, ty, len)

vec = [tx(ctr(1),ctr(2)),ty(ctr(1),ctr(2))];

cur_pt = ctr;
x = ctr(2);
y = ctr(1);

n = 1;
while (n<=len)
    cur_pt = round(cur_pt+[ty(cur_pt(1),cur_pt(2)),tx(cur_pt(1),cur_pt(2))]);
    n = n+1;
    if cur_pt(1)<1 || cur_pt(2)<1 || cur_pt(1)>size(tx,1) || cur_pt(2)>size(tx,2) 
        break; % outside the boundary of the image
    end 
    if (tx(cur_pt(1),cur_pt(2)) == 0 && ty(cur_pt(1),cur_pt(2))==0)
        break; % direcitonal vector equals to zero
    end
    if sum(x==cur_pt(2) & y==cur_pt(1))
        break; % same point is found which means a loop is formed
    end
    x = [x; cur_pt(2)];
    y = [y; cur_pt(1)];
end
n = 1;
cur_pt = ctr;
while (n<=len)
    cur_pt = round(cur_pt-[ty(cur_pt(1),cur_pt(2)),tx(cur_pt(1),cur_pt(2))]);
    n = n+1;
    if cur_pt(1)<1 || cur_pt(2)<1 || cur_pt(1)>size(tx,1) || cur_pt(2)>size(tx,2) 
        break; % outside the boundary of the image
    end 
    if (tx(cur_pt(1),cur_pt(2)) == 0 && ty(cur_pt(1),cur_pt(2))==0)
        break; % direcitonal vector equals to zero
    end
    if sum(x==cur_pt(2) & y==cur_pt(1))
        break; % same point is found which means a loop is formed
    end
    x = [cur_pt(2); x];
    y = [cur_pt(1); y];
end

[~,pos] = min((x-ctr(2)).^2+(y-ctr(1)).^2);