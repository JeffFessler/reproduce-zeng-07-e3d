function [cost, xgrad,ygrad, zgrad] = calc_corr_grad(ob_bsp, theta_x, theta_y, theta_z,...
                  g, g_bar,img_coeff, xxi,yyi,zzi, Gb, isgrad)

if(nargin==11) isgrad=1; end;


parasize=length(theta_x);
na=size(g,3);
[xsize,ysize,zsize]=size(img_coeff);
vol=xsize*ysize*zsize;
Bx=ob_bsp.Bx; By=ob_bsp.By;Bz=ob_bsp.Bz;Bt=ob_bsp.Bt;
xgrad=zeros(size(theta_x)); ygrad=xgrad; zgrad=xgrad;

d1=0;
for i=1:na

geomx = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_x));
geomy = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_y));
geomz = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_z));

%image interpolation
img_out=single(zeros(vol,1));
imggradx=img_out; imggrady=img_out; imggradz=img_out;

%[xxi,yyi,zzi]=ndgrid(xi,yi,zi);
    [img_out,imggradx, imggrady, imggradz]=...
        interpolate(img_coeff,double(reshape(geomx+xxi,vol,1)),...
        double(reshape(geomy+yyi,vol,1)),...
        double(reshape(geomz+zzi,vol,1)), 3, 1,1,1);

%clear xxi yyi zzi;
img_out=single(reshape(img_out,xsize,ysize,zsize));
img_out(find(img_out<0))=0;

imggradx=single(reshape(imggradx,xsize,ysize,zsize));
imggrady=single(reshape(imggrady,xsize,ysize,zsize));
imggradz=single(reshape(imggradz,xsize,ysize,zsize));

%bp_gp = A'((p_m - p_bar)/T1 - (p_m - p_bar)/T2)
p=Gb{i}*img_out;
d1=d1 - log(corr2(g(:,:,i),p));
if(isgrad)
    p_bar(i) = mean(p(:));
    T1 = sum(sum((g(:,:,i)-g_bar(i)).*(p-p_bar(i))));
    T2(i) = sum(sum((p-p_bar(i)).^2));  
    bp_gp = Gb{i}'*((p-p_bar(i))/T2(i)-(g(:,:,i)-g_bar(i))/T1);
%gradient
    xgrad=reshape(kron_product4(Bx', By', Bz', Bt(i,:)', (bp_gp).*imggradx ), parasize,1) + xgrad;
    ygrad=reshape(kron_product4(Bx', By', Bz', Bt(i,:)', (bp_gp).*imggrady ), parasize,1) + ygrad;
    zgrad=reshape(kron_product4(Bx', By', Bz', Bt(i,:)', (bp_gp).*imggradz ), parasize,1) + zgrad;
    clear bp_gp;
end;


clear imggradx imggrady imggradz img_out;
end;
cost=d1/na;
xgrad=xgrad/na;
ygrad=ygrad/na;
zgrad=zgrad/na;
