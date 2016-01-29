function b = Norm(a,mask,para)

if nargin < 3
    para = 'default';
    if nargin < 2
        mask = ones(size(a));
    end
end

switch lower(para)
    case 'zeromean'
        mu = mean(a(mask==1));
        sigma = std(a(mask==1));
        b = (a-mu)./sigma;
    otherwise
        b = (a-min(a(mask==1)))./(max(a(mask==1))-min(a(mask==1)));
end
