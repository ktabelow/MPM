% 
% function to smooth a three-dimensional array given a bandwidth
%
% =========================================================================
%
% 2017/03/24
%
% by C. D'Alonzo
% matlab version of the ckernsm function in R package aws
%
% The function smooth3d_with_kern(y, h) calculates a three-dimensional
% smoothing of the 3D array y using the bandwidth h and an epanechnicov
% kernel (supported on (-1,1)) on a grid
%
% Input:
%   y          - a 3D array to be smoothed
%   h          - bandwidth for the smoothing
% Output:
%   y_smoothed - the 3D smoothed version of y
function [y_smoothed] = smooth3d_with_kern(y, h)

    
        
dy = size(y);
d = length(dy);
if d~=3
    error('y has to be a 3D array, but it has size %d',d);
end

if length(h)==1
    h = ones(1,d).*h;
end

if length(h)~=d
    error('Incompatible length of bandwidth vector h.\nHas to be a scalar or a vector of length 3, but it was %d',length(h));
end

y = reshape(y,size(y,1),size(y,2)*size(y,3));
 
z = extend_y(y,h(1));
yy = z.yy;
dyy1 = size(yy,1);

kern1 = max(0, 1-grid_for_smoothing(dyy1,h(1)).^2)/(h(1)*4/3);
f1 = fft(yy);
f2 = fft(kern1);
invFourier = ifft(f1.*(f2'*ones(1,size(f1,2))));
yhat = invFourier(z.ind,:)/sum(kern1);
yhat = reshape(yhat,dy);

yhat = permute(yhat,[2 1 3]);

yhat = reshape(yhat,dy(2),dy(1)*dy(3));


z = extend_y(yhat,h(2));
yy = z.yy;
dyy2 = size(yy,1);


kern2 = max(0, 1-grid_for_smoothing(dyy2,h(2)).^2)/(h(2)*4/3);
f1 = fft(yy);
f2 = fft(kern2);
invFourier = ifft(f1.*(f2'*ones(1,size(f1,2))));
yhat = invFourier(z.ind,:)/sum(kern2);

yhat = reshape(yhat, dy(2),dy(1),dy(3));

yhat = permute(yhat,[3 2 1]);

yhat = reshape(yhat, dy(3),dy(1)*dy(2));

z = extend_y(yhat,h(3));
yy = z.yy;
dyy3 = size(yy,1);
kern3 = max(0, 1-grid_for_smoothing(dyy3,h(3)).^2)/(h(3)*4/3);
f1 = fft(yy);
f2 = fft(kern3);
invFourier = ifft(f1.*(f2'*ones(1,size(f1,2))));
yhat = invFourier(z.ind,:)/sum(kern3);

yhat = reshape(yhat, dy(3),dy(1),dy(2));
yhat = permute(yhat,[2 3 1]);

y_smoothed = yhat;




function [gd] = grid_for_smoothing(d,h)

d0 = floor(d/2)+1;
%step = 1/(d0-1);
%gd = 0:step:1;
gd = linspace(0,1, d0);
if 2*d0 == d+1
    gd = [gd -gd(d0:-1:2)];
else
    gd = [gd -gd((d0-1):-1:2)];
end

gd = gd/2/h*d;


function [z] = extend_y(y,h)

n = size(y,1);
h = min(h, floor(n/2));
nn = nextn(n+2*h);
yy = zeros(size(y,2),nn);
ih0 = floor((nn-n)/2);
ih1 = nn-ih0-n;
ind = (ih0+1):(ih0+n);
yy(:,ind) = y';
yy(:,1:ih0) = y(ih0:-1:1,:)';
yy(:,(nn-ih1+1):nn) = y(n:-1:(n-ih1+1),:)';
z.yy = yy';
z.ind = ind;

function [nextComposite] = nextn(n)
nextComposite = floor(n);

while nextComposite<10000000 
    
    result = nextComposite;
    while mod(result,2)==0 
        result = floor(result/2);        
    end
    while mod(result,3)==0 
        result = floor(result/3);        
    end
    while mod(result,5)==0 
        result = floor(result/5);        
    end
    if result==1 
        break;
    end
    nextComposite = nextComposite +1;
end



    





