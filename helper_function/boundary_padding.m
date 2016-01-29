function [im_out, mask_out] = boundary_padding(im_in, im_mask, choice, offset)

if nargin < 3
    choice = 'simple';
end

if nargin < 4
    offset = 0;
end

disp('Start padding...');
tic
switch lower(choice)
    case 'simple'
        %% boundary padding
        dim = size(im_in);
        N = round(min(dim)/10);
        psz = round(max(dim)/10);
        [rm, cm] = find(im_mask);
        idm = find(rm==cm, 1);
        sub_im = im_in(rm(idm):rm(idm)+psz,cm(idm):cm(idm)+psz);

        im_temp = repmat(sub_im, ceil((dim(1)+N*2)/psz), ceil((dim(2)+N*2)/psz));
        im_temp = im_temp(1:dim(1)+N*2,1:dim(2)+N*2);

        im_large = zeros(dim+N*2);
        im_large(N+1:end-N,N+1:end-N) = im_mask;
        lgc = (im_large==0);
        im_large(N+1:end-N,N+1:end-N) = im_in;
        im_large(lgc) = im_temp(lgc);
        
    case 'draw'
        dim = size(im_in);
        N = 20;
        
        se = strel('disk',1);
        im_large = zeros(dim+N*2);
        im_large_mask = zeros(dim+N*2);
        im_large_mask(N+1:end-N,N+1:end-N) = im_mask;     
        if offset
            im_large_mask = imerode(im_large_mask,se);
            im_large_mask = imerode(im_large_mask,se);
            im_large_mask = imerode(im_large_mask,se);
        end
        bw_bd = im_large_mask-imerode(im_large_mask,se);
        im_large(N+1:end-N,N+1:end-N) = im_in;
      
        n = N;
        
        [r,c] = find(bw_bd);
        queue = [r,c];
        q_idx = find(bw_bd);
        while n
            outer_bdy = imdilate(im_large_mask,se) - im_large_mask;
            im_large_mask = imdilate(im_large_mask,se);
            [or, oc] = find(outer_bdy==1);
            [~, idx] = min(sqrt(bsxfun(@minus, or, queue(:,1)').^2+bsxfun(@minus, oc, queue(:,2)').^2),[],2);
            im_large(outer_bdy==1) = im_large(q_idx(idx));
            queue = [or, oc];
            q_idx = sub2ind(size(im_large),queue(:,1),queue(:,2));
            n = n-1;
        end
        
    otherwise
        error('Unidentified choice!');
end

im_out = im_large;
mask_out = zeros(size(im_large));
mask_out(N+1:end-N,N+1:end-N) = 1;

disp('...Finish padding');
toc