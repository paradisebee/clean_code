function pairwise = setSmoothness2Dfunc(im,lumbda,sigma,choice,mask,func_handle)

switch choice
    case '4nbh'
        flag = 1;
    case '8nbh'
        flag = 2;
end

if nargin < 6
    func_handle = @(a,b) exp(-(a-b)^2/2/sigma^2);
end
if nargin < 5
    mask = ones(size(im));
end


[row,col] = size(im);

L = row*col;

indexArray1 = [];
indexArray2 = [];
valueArray = [];
for m = 1:L
    % 4 neighborhood
    if m-row > 0,              % i, j-1
        indexArray1(end+1) = m;
        indexArray2(end+1) = m-row;
        if mask(m-row)==0 || mask(m) == 0
            valueArray(end+1) = 1;
        else
            valueArray(end+1) = func_handle(im(m),im(m-row));
        end
    end
    
    if mod(m-1,row)~=0,        % i-1, j        
        indexArray1(end+1) = m;
        indexArray2(end+1) = m-1;
        if mask(m-1)==0 || mask(m) == 0
            valueArray(end+1) = 1;
        else
            valueArray(end+1) = func_handle(im(m),im(m-1));
        end
    end
    
    if m+row <= L,             % i, j+1
        indexArray1(end+1) = m;
        indexArray2(end+1) = m+row;
        if mask(m+row)==0 || mask(m) == 0
            valueArray(end+1) = 1;
        else
            valueArray(end+1) = func_handle(im(m),im(m+row));
        end
    end
    
    if mod(m,row)~=0,          % i+1, j
        indexArray1(end+1) = m;
        indexArray2(end+1) = m+1;
        if mask(m+1)==0 || mask(m) == 0
            valueArray(end+1) = 1;
        else
            valueArray(end+1) = func_handle(im(m),im(m+1));
        end
    end

    % 8 neighborhood
    if flag == 2
        if m-row-1 > 0 && mod(m-row-1,row)~=0, % i-1, j-1
            indexArray1(end+1) = m;
            indexArray2(end+1) = m-row-1;
            if  mask(m-row-1)==0 || mask(m) == 0
                valueArray(end+1) = 1;
            else
                valueArray(end+1) = func_handle(im(m),im(m-row-1));
%             valueArray(end+1) = 1/sqrt(2)*func_handle(im(m),im(m-row-1));
            end
        end
        if m-row+1 > 0 && mod(m-row,row)~=0,  % i+1, j-1
            indexArray1(end+1) = m;
            indexArray2(end+1) = m-row+1;
            if mask(m-row+1)==0 || mask(m) == 0
                valueArray(end+1) = 1;
            else
                valueArray(end+1) = func_handle(im(m),im(m-row+1));
%             valueArray(end+1) = 1/sqrt(2)*func_handle(im(m),im(m-row+1));
            end
        end
        if m+row-1 <= L && mod(m+row-1,row)~=0, % i-1, j+1
            indexArray1(end+1) = m;
            indexArray2(end+1) = m+row-1;
            if mask(m+row-1)==0 || mask(m) == 0
                valueArray(end+1) = 1;
            else
                valueArray(end+1) = func_handle(im(m),im(m+row-1));
%             valueArray(end+1) = 1/sqrt(2)*func_handle(im(m),im(m+row-1));
            end
        end
        if m+row+1 <= L && mod(m+row,row)~=0, % i+1, j+1
            indexArray1(end+1) = m;
            indexArray2(end+1) = m+row+1;
            if mask(m+row+1)==0 || mask(m) == 0
                valueArray(end+1) = 1;
            else
                valueArray(end+1) = func_handle(im(m),im(m+row+1));
%             valueArray(end+1) = 1/sqrt(2)*func_handle(im(m),im(m+row+1));
            end
        end
    end
        
end
valueArray = lumbda * valueArray;
pairwise = sparse(indexArray1,indexArray2,valueArray,row*col,row*col);

