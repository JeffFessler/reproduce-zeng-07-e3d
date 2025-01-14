function [mask,I,J,K]=find_mask(img,level);
%[mask,I,J,K]=find_mask(img, level);

mask=zeros(size(img));
mask(find(img>level))=1;
mask=logical(mask);
for i=1:size(mask,3)
   mask(:,:,i)=imfill(mask(:,:,i),'holes');
end;

ind=find(mask==1);
dim=length(size(img));
if(dim==2)
    [I,J]=ind2sub(size(img),ind);
    K=0;
end
if(dim==3)
    [I,J,K]=ind2sub(size(img),ind);
end
