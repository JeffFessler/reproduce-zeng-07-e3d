function [imgout, mask] = rm_air(imgin,N)
%[imgout, mask] = rm_air(imgin,N)
%N=2,3,4 etc

[nx,ny,nz]=size(imgin);
x = grayslice(imgin./max(imgin(:)),N);
x(find(x>0))=1;
for i=1:nz
    x(:,:,i) = imfill(x(:,:,i),'holes');
end;
imgout = imgin.*double(x);
%[xc,yc] = find(x==1);
%xrange = [min(xc) max(xc)];
%yrange = [min(yc) max(yc)];
mask = logical(x);