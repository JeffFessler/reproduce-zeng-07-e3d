function hess = calc_cor_hess_a(ob_bsp, theta_x, theta_y, theta_z,...
                  img_coeff, na, xxi,yyi,zzi, Gb, d)

parasize=length(theta_x);

[xsize,ysize,zsize]=size(img_coeff);
vol=xsize*ysize*zsize;
Bx=ob_bsp.Bx; By=ob_bsp.By;Bz=ob_bsp.Bz;Bt=ob_bsp.Bt;
xgrad=zeros(size(theta_x)); ygrad=xgrad; zgrad=xgrad;

hess=0;
for i=1:na
voldir = zeros(xsize,ysize,zsize);

geomx = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_x));
geomy = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_y));
geomz = kron_product4_1d44(single(Bx), single(By), single(Bz),single(Bt(i,:)),single(theta_z));

%image interpolation
img_out=single(zeros(xsize*ysize*zsize,1));
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
%
p=Gb{i}*img_out;
nb=length(p(:));
p_bar=mean(p(:));
T=sum((p(:)-p_bar).^2)/2;

dirx = kron_product4_1d44(Bx, By, Bz, Bt(i,:), d(1:parasize));
       %single(kron_product4_d4(xsize, ysize, zsize, Bt(i,:), partial_dirx));
voldir=single(imggradx.*dirx); clear imggradx dirx;
diry = kron_product4_1d44(Bx, By, Bz, Bt(i,:), d([1:parasize]+parasize));
       %single(kron_product4_d4( xsize, ysize, zsize, Bt(i,:), partial_diry));
voldir=single(voldir+imggrady.*diry); clear imggrady diry;
dirz = kron_product4_1d44(Bx, By, Bz, Bt(i,:), d([1:parasize]+2*parasize));
       %single(kron_product4_d4(xsize, ysize, zsize, Bt(i,:), partial_dirz));
voldir = single(voldir+imggradz.*dirz); clear imggradz dirz;
 
proj_voldir = calc_fp(Gb,voldir,i); clear voldir;

hess = (sum(double(proj_voldir(:)).^2)-sum(double(proj_voldir(:)))^2/nb)/T + hess;

end;
hess=hess/na;