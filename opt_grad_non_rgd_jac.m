%function [final_jacmin, final_cor]=coropt_grad_non_rgd_jac(gain,reg_param) 
clear all;
test=0
path(path,'../');
%step=2
step=input('Multi-resolution Step:');

file1=input('File 1 (reference):\n','s'); 
file2=input('File 2 (Target):\n','s');
tmpfile=input('Temporary folder to store the coarser reg results:\n','s');
ctex=fld_read(file1);
ctin=fld_read(file2);
ctex=downsample3d(ctex,2,1);
ctin=downsample3d(ctin,2,1);

%downsample of the images
ds_x = 2/2^(step-1); %ds_y = 8; ds_z = 8;
gain=10^5*5^(step-1)%5*10^6*(5^(step-1)) for phantom    


sub_smp_xax=8;%16;%4;%16;%8;%1;
sub_smp_yax=8;%16;%4;%16;%8%1;
sub_smp_zax=4;%8;%4;%8;%6;%1


new_ct1_d4=downsample3(ctex,ds_x);%ct20_lp(1:ds_x:end,1:ds_y:end,1:ds_z:end);
new_ct2_d4=downsample3(ctin,ds_x);%ct80_lp(1:ds_x:end,1:ds_y:end,1:ds_z:end);

mask=find_mask(new_ct2_d4,5);
clear ctex ctin ct20 ct20_lp ct80 ct80_lp;

coeff=coeff_3d_mex(single(new_ct1_d4),3);
coeff=double(coeff);

[x_size y_size z_size]=size(coeff);


%--------------------------------------------------
%  Preparing sub-image here 
%--------------------------------------------------
image=double(new_ct2_d4);

x_bnd=0; y_bnd=0; z_bnd=0;

xi=linspace(0, x_size-1, x_size);
yi=linspace(0, y_size-1, y_size);
zi=linspace(0, z_size-1, z_size);

ref=image;


%------------------------------------
% Setting B-spline function Sets 
%-----------------------------------
x_st=1;
[aux x_ed]=size(xi);

y_st=1;
[aux y_ed]=size(yi);

z_st=1;
[aux z_ed]=size(zi);

if step==1
    x_locx=[x_st-1:sub_smp_xax:x_ed];
    x_locy=[y_st-1:sub_smp_yax:y_ed];
    %x_locz=[z_st-1:sub_smp_zax:z_ed];
    x_locz= [z_ed:-sub_smp_zax:z_st-1];

    y_locx=x_locx;
    y_locy=x_locy;
    y_locz=x_locz;

    z_locx=x_locx;
    z_locy=x_locy;
    z_locz=x_locz;
    [y_x y_y y_z]=meshgrid(y_locx, y_locy,y_locz);

size_parm=length(y_locx)*length(y_locy)*length(y_locz);
xcoeff=zeros(size_parm,1);
ycoeff=zeros(size_parm,1);
zcoeff=zeros(size_parm,1);
end

%-----------------------------------------
%initialize bspline coefficients
%-----------------------------------------

if(step>1)
    inifile=sprintf('temp%s_new',num2str(step-1));
    load([tmpfile '/' inifile]);
    y_locx=new_locx; y_locy=new_locy; y_locz=new_locz;
    sub_smp_xax = y_locx(2)-y_locx(1);
    sub_smp_xay = y_locy(2)-y_locy(1);
    sub_smp_xaz = y_locz(2)-y_locz(1);
    [y_x y_y y_z]=meshgrid(y_locx, y_locy,y_locz);
    size_parm=length(y_locx)*length(y_locy)*length(y_locz);
end;


%-----------------------------------------------------------------------------
[y_x y_y y_z]=meshgrid(y_locx, y_locy,y_locz);
i=0;

lambda=0.000;
hesst=[];
gradt=[];

theta_x=xcoeff(:);
theta_y=ycoeff(:);
theta_z=zcoeff(:);



%gain = alpha
reg_param = 8*1e5;
slop = 1e2;

diff_cost= 10^-3;
diff_cort=51;

jac_thr=0.05;

cort=0;
jac_min=1;
lev = [4*0.1882 4*0.1882 1*0.5]; %voxel resolution
max_dif = 100;
i=1;
% load ~rzeng/MotionImg/phan/results/reg_ct02_ds221_16168.mat;
%load ~rzeng/4dct/result/4287/reg1.mat;
tic
while (jac_min<0 | diff_cost>5*10^(-7)/ds_x) %max_dif>0.1);%min_dif>0.1)%cort<0.96)% 

if(step>1 & i>30) break; end
if(step>1 & i>60) break; end

[f, img, geomx, geomy, geomz, pen,  grx, gry, grz, hess_x, hess_y, hess_z, jacobian] = feval('fob_non_rgd_jac',theta_x, theta_y, theta_z, y_x, y_y, y_z, ref, coeff, lambda, reg_param, slop,sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, xi, yi, zi, jac_thr,lev);

jac_min=min(min(min(jacobian)))

if i>1 jacmin(i)=jac_min; end;

%-------------------------------------------------------------------------------
%
% Update routine by gradient descent
%-------------------------------------------------------------------------------

ix=find(hess_x);
iy=find(hess_y);
iz=find(hess_z);
theta_old = [theta_x;theta_y;theta_z];

theta_x = theta_x  - gain * grx; %  ./ hess_x(ix);
theta_y = theta_y  - gain * gry; %  ./ hess_y(iy);
theta_z = theta_z  - gain * grz; %  ./ hess_z(iz);
if (i>1)
 theta_new = [theta_x;theta_y; theta_z];
 theta_dif = theta_new - theta_old;
 max_dif = max(abs(theta_dif(:)));%/mean(abs(theta_old))*100;
 maxdif(i) = max_dif;
end;    

    fvalt(i)=f;


    if (i > 1)
       diff_cost(i) = fvalt(i-1) - fvalt(i)    
    %diff_per(i) =  (fvalt(i-1) - fvalt(i))/
       if(diff_cost(i)<0) 
           disp 'unstable, please decrease the value of "gain".'; 
           break ; 
       end
    end

    
    corr(i)=corr2(img(:), ref(:));%corr2(img(mask(:)), ref(mask(:)))
    cort = corr(i)
    pentt(i)=pen;
    jactt(i)=jac_min;    
    
i=i+1
%    time=cputime-time
end
toc
%save temp1.mat theta_x theta_y theta_z x_locx x_locy x_locz corr...
% fvalt img jacobian jacmin difmin;
theta_x=theta_old(1:size_parm); theta_y=theta_old([1:size_parm]+size_parm);
theta_z=theta_old([1:size_parm]+2*size_parm);
%save temp.mat theta_x theta_y theta_z geomx geomy geomz y_locx y_locy y_locz;
%save geo.mat geomx geomy geomz;

beep;
if(step<2) run b_expand; end;
