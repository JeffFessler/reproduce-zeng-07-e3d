%-----------------------------------------------------
%Read the reference images and store them in one file
%
%Inputs
%  str1: the directory that contains all the CT slices
%  mu_h2o: the attenuation coefficient of water, can be changed by the user depending 
%        on the X-ray energy.
%  str2: the location to store the reference CT volume
%
% User specified
%   The table romoving and image cropping may need to be modified based on the images.
%-----------------------------------------------------

str1='/n/ir24/y/rzeng/dov/newDOVexperiment/staticCT/pos1/';
str2='/n/ir24/y/rzeng/dov/temp/phan/ct1.fld';
mu_h2o = 0.0158; 

D=dir(str1);
nx=512; ny=512; nz=length(D)-2;
x=zeros(nx, ny, nz);
for i=1:nz
    str=[str1 D(i+2).name];
    img=dicomread(str);
    x(:,:,i)=img;
end;
x(find(x==-2000))=0;
x=HU2mu(x,mu_h2o);
x=single(x);

%remove table and crop the image--- Need to modify it base on the images
tblmask=ones(size(x));
t=[1:512]; s=((t-256)/256).^2*20;
sta=512-120-round(s);
for i=1:512 tblmask(sta(i):end,:,:)=0; end;
tblmask=logical(tblmask);
y=x.*tblmask;
ct=single(y(111:end-110, 91:end-90,:));

fld_write(str2,ct);
clear img y z;