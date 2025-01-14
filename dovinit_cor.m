%-----------------------------------------------
%DOV with the proportionality motion model: estimating alpha_t[1xna], t=1:na,
%
%Inputs 
%   ctfile:reference volume.
%   Projfile:projection views, 
%   geomx, geomy, geomz: extreme deformation fields
%
%output
%   alpha_t[1xna]: (saved in file 'a.mat')
%------------------------------------------------
clear all;
step=1;

%The reference volume
ctfile='/n/ir24/y/rzeng/dov/phan/ct2_aff3.fld';
ct=fld_read(ctfile);

%The cone-beam projection views
projfile='/n/ir24/y/rzeng/dov/phan/dynali_aperiod.fld'
proj=fld_read(projfile);%proj_deformedct2_sim.fld
g=proj(end:-1:1,:,:);%(:,:,1:2:end); 
%g=g(end:-1:1,:,end:-1:1); %for ct0 CBCT only
clear proj;

%downsampling
if(1)
dsx=4/(2^(step-1));%dsy=2;dsz=2;
dsz=2/(2^(step-1));
ct=downsample3d(ct,dsx,dsz);%d(ct,dsx,dsz);
ct=single(ct);%*10;

dsp=2/(2^(step-1));
g=downsample3d(g,dsp,1);%/dsx;

 g=single(g);
[nu,nv,na]=size(g);

end;

%Image interpolation
img_coeff=coeff_3d_mex(single(ct),3);
img_coeff=double(img_coeff);


[nx,ny,nz]=size(ct);
na=size(g,3);
mot_num = na;


[xsize ysize zsize]=size(ct);

xi=linspace(0,xsize-1,xsize);
yi=linspace(0,ysize-1,ysize);
zi=linspace(0,zsize-1,zsize);
ti=linspace(0,mot_num-1, mot_num);  


% The full deformation
load data/geo_ct02_16168.mat;
geomx=downsample3d(geomx,dsx,dsz)/dsx;geomy=downsample3d(geomy,dsx,dsz)/dsx;geomz=downsample3d(geomz,dsx,dsz)/dsz;
% geomx=-geomx;geomy=-geomy;geomz=-geomz;

%---smoothness penalty----------------
ele=[1 -2 1 zeros(1,na-3)];  %For 2nd order penalty
C1a = zeros(na);
for i=2:na-1
    C1a(i,:)=wshift('1D',ele,-(i-2));
end;
C1a(1,:)=[-1 1 zeros(1,na-2)];
C1a(end,:)=[zeros(1,na-2) 1 -1];
C2a=C1a'*C1a;

%--- For projection--------
run cone;%
g_bar=mean(reshape(g, nu*nv, na));
%-----------------------------------

[xxi,yyi,zzi]=ndgrid(xi,yi,zi);

vol=xsize*ysize*zsize;


%---optimization--------------------
tolf=1e-5;%1e-5;%1e-6;%10;
tolch=0.01;
beta= 0.01;% The roughness penalty parameter
thr = 7;
fold=1000;
fdif=tolf+1;
max_deltaa=0.1;

iter=1;
d1=3000;
reset = 0;
a=zeros(na,1);
H=zeros(na); 
%load a.mat; a=a*2;
while (iter<30 || max_deltaa>tolch)%(abs(fdif)>tolf)%(d1>2100)%

    iter

 
d1=0; 
tic


for i=1:na

%image interpolation
img_out=single(zeros(xsize*ysize*zsize,1));
imggradx=img_out; imggrady=img_out; imggradz=img_out;
 
%[xxi,yyi,zzi]=ndgrid(xi,yi,zi);
    [img_out,imggradx, imggrady, imggradz]=...
        interpolate(img_coeff,double(reshape(a(i)*geomx+xxi,vol,1)),...
        double(reshape(a(i)*geomy+yyi,vol,1)),...
        double(reshape(a(i)*geomz+zzi,vol,1)), 3, 1,1,1);

%clear xxi yyi zzi;
img_out=single(reshape(img_out,xsize,ysize,zsize));
%img_out(find(img_out<0))=0;
imggradx=single(reshape(imggradx,xsize,ysize,zsize));
imggrady=single(reshape(imggrady,xsize,ysize,zsize));
imggradz=single(reshape(imggradz,xsize,ysize,zsize));

%Calculating the cost function and the derivatives
p=Gb{i}*img_out;
d1=d1 - log(corr2(g(:,:,i),p));

p_bar(i) = mean(p(:));
T1 = sum(sum((g(:,:,i)-g_bar(i)).*(p-p_bar(i))));
T2(i) = sum(sum((p-p_bar(i)).^2));  
gp = ((p-p_bar(i))/T2(i)-(g(:,:,i)-g_bar(i))/T1);
%(Af-P)
temp1=imggradx.*geomx+imggrady.*geomy+imggradz.*geomz;
temp2=Gb{i}*temp1;
agrad(i)=gp(:)'*temp2(:);
H(i,i) = (temp2(:)'*temp2(:)-sum(temp2(:))^2/nv/nu)/T2(i);
clear pe;

clear imggradx imggrady imggradz img_out;
end;
agrad=agrad/na;
H=H/na;
%ahess=ahess/nu/nv/na;
%bhess=bhess/nu/nv/na;
%chess=chess/nu/nv/na;
%     % Mean Square error
    disp('The cost is :..............');
    d1=d1/na %mean(mean((ref_sino-sino_dyna).^2));
    dd(iter)=d1;
   
    %------penalty-----------------------------
    apen = (C1a*a)'*(C1a*a)/na;
    apen_grad = (C1a'*C1a)*a/na;
    pen_hess = (C1a'*C1a)/na;
    pen = beta*(apen);
    
    dist(iter)=d1+pen;
    dist(iter)
    fdif(iter)=dist(iter)-fold;
   % fold=dist(iter); 
    if(fdif(iter)>0) disp 'reset'; reset=1; fold=dist(iter); continue; end;
    
     grad = [agrad(:)]+ [apen_grad(:)]*beta;
 %    hess = [ahess(:); bhess(:); chess(:)]+ [pen_hess(:);pen_hess(:);pen_hess(:)]*beta;
     theta = [a;];
     
    %Conjugate gradient update 
     p = grad; %pp*grad;
     if (iter==1|reset==1)
         r=0; d_old=zeros(na,1);
         reset=0;
     else
         r=real(grad'*(grad -grad_old)) / real(grad_old'*grad_old);
     end
     d = -grad + r*d_old;             
         
     theta_old=theta;
     alpha0=0;
  
     alpha_grad = d'*grad;
     hess=H+beta*pen_hess;
     alpha_hess = d'*hess*d;
     alpha0 = -alpha_grad/alpha_hess
   
   if(alpha0<0) %reset cg
     d_old=0;
     reset=1;
     fold=dist(iter); 
     continue; %fold=dist(iter-1);
   else
   
        dir = alpha0*d;     
        theta = theta_old + dir;           
        a = theta(1:na);

        grad_old=grad;
        
        d_old = d;
        max_deltaa=max(abs(dir))
   end;
   fold=dist(iter); 
    if (iter>50) break; end;  
    
    iter=iter+1;
     
     beep
     toc
end;
save a.mat a;

beep;


clear partial_dirx partial_diry partial_dirz partial_geomx partial_geomy partial_geomz;

% get img_out back;
% [xxi,yyi,zzi]=ndgrid(xi,yi,zi);
% for i=1:na
%    
%     [img_est(:,i),imggradx, imggrady, imggradz]=...
%         interpolate(img_coeff,double(reshape(a(i)*geomx+xxi,vol,1)),...
%         double(reshape(a(i)*geomy+yyi,vol,1)),...
%         double(reshape(a(i)*geomz+zzi,vol,1)), 3, 1,1,1);
% 
% end;
% clear xxi yyi zzi imggradx imggrady imggradz;
% img_est=single(reshape(img_est,nx,ny,nz,mot_num));
% img_est(find(img_est<0))=0;
