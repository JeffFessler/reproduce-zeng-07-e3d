function [i,j,k,l]=find_sub(x,val)
%[i,j,k,l]=find_sub(x,val)

id=find(x==val);
[i,j,k,l]=ind2sub(size(x),id);
