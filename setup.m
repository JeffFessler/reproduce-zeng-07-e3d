
%step=2;

ctex=fld_read('/n/ir24/y/rzeng/dov/phan/ct2_aff3.fld');
proj=fld_read('/n/ir24/y/rzeng/dov/phan/dynali_aperiod.fld');
proj_dyna=proj(end:-1:1,:,1:80); 
%proj_dyna=proj_dyna(end:-1:1,:,:);%(end:-1:1,:,end:-1:1)for ct0 CBCT only
clear proj;


if step==1 
%cyc=[1 16; 20 36; 40 56; 60 82];%Oracle (best)
%cyc=[0 16; 21 33; 37 48; 55 69]; %manual
%cyc=[3 18; 22 36; 37 52; 57 69]; %corr 
%auto placement
cyc= find_breathcyc
end
knotpercyc=5;

%downsampling
if(1)
    dsx=4/(2^(step-1));%dsy=2;dsz=2;
    dsz=2/(2^(step-1));
    ct=downsample3d(ctex,dsx,dsz);
    ct=single(ct); clear ctex;

    dsp=4/(2^(step-1));
    g=downsample3d(proj_dyna,dsp,1);%/dsx;
%proj_dyna=reshape(proj_dyna,size(proj_dyna,1)*size(proj_dyna,2),size(proj_dyna,3)); %for 3l
    g=single(g); clear proj_dyna;
%crop the rigid plane
    %nv=size(g,2);
    %g=g(:,1:floor(nv*4/5),:);
end;
[nu,nv,na]=size(g);


img_coeff=coeff_3d_mex(single(ct),3);
img_coeff=double(img_coeff);


sub_smpx =16;%16;%32;%16;%8;%4;
sub_smpy =16;%16;% 32;%16;%8;%4;
sub_smpz =8;%8;%16;%8;%4;
%sub_smpt = 4;%2;%4;
[nx,ny,nz]=size(ct);
na=size(g,3);
mot_num = na;


[xsize ysize zsize]=size(ct);

xi=linspace(0,xsize-1,xsize);
yi=linspace(0,ysize-1,ysize);
zi=linspace(0,zsize-1,zsize);
ti=linspace(0,mot_num-1, mot_num);  

mask = find_mask(ct,0);

  run placeknots;
  
[xlocx xlocy xlocz xloct]=ndgrid(x_x,x_y, x_z,x_t); %xlocz = 0;

ylocx=xlocx; ylocy=xlocy; ylocz=xlocz; yloct=xloct;
zlocx=xlocx; zlocy=xlocy; zlocz=xlocz; zloct=xloct;

parasize = length(x_x)*length(x_y)*length(x_z)*length(x_t)
theta_x = zeros(1,parasize)';
theta_y = zeros(1,parasize)';
theta_z = zeros(1,parasize)';


%---smoothness penalty----------------
clt = length(x_t);  lx = length(x_x); ly = length(x_y); lz = length(x_z);
ltt=length(x_t);
%ele=[-1 1 zeros(1,mot_num-2)];   % For  1st order penalty   
ele=[1 -2 1 zeros(1,ltt-3)];  %For 2nd order penalty
C1t = zeros(ltt);
for i=2:length(x_t)-1
    C1t(i,:)=wshift('1D',ele,-(i-2));
end;
ele=[1 -2 1 zeros(1,lx-3)];
C1x = zeros(lx);
for i=2:length(x_x)-1
    C1x(i,:)=wshift('1D',ele,-(i-2));
end;
ele=[1 -2 1 zeros(1,ly-3)];
C1y = zeros(ly);
for i=2:length(x_y)-1
    C1y(i,:)=wshift('1D',ele,-(i-2));
end;
ele=[1 -2 1 zeros(1,lz-3)];
C1z = zeros(lz);
for i=2:length(x_z)-1
    C1z(i,:)=wshift('1D',ele,-(i-2));
end;
C2t=C1t'*C1t;
C2x=C1x'*C1x;
C2y=C1y'*C1y;
C2z=C1z'*C1z;

%temporal aperiodicity penalty
full=1;
%full aperiodicity penalty
if (full)
  %elep = [-1 zeros(1, (lt-1)/2-1) 1 zeros(1,(lt-1)/2)];%exact T  
  elep = [-1 zeros(1, knotpercyc-1) 1 zeros(1,ltt-knotpercyc-1)];
  %C1p = zeros((lt-1)/2+1, lt);
  for i=1:knotpercyc*(length(x_t)/knotpercyc-1) %1:(lt-1)/2+1
      C1p(i,:) = wshift('1D',elep, -(i-1)); 
  end
  C1px=C1p;C1py=C1p;C1pz=C1p;
%Penalize peaks only
else
  elep = [0 0 1 0 0  0 0 -1 0 0 zeros(1,knotpercyc*(size(cyc,1)-2))];
  %C1p = zeros((lt-1)/2+1, lt);
  for i=1:(size(cyc,1)-1) %1:(lt-1)/2+1
      C1p(i,:) = wshift('1D',elep, -(i-1)*knotpercyc); 
  end
  C1px=[0 0 1/2 0 0  0 0 -1 0 0  0 0 0 0 0  0 0 1/2 0 0;
       0 0 1/2 0 0  0 0 0 0 0   0 0 -1 0 0 0 0 1/2 0 0];
  C1py=[0 0 -1 0 0  0 0 1/2 0 0  0 0 1/2 0 0  0 0 0 0 0;
        0 0 0 0 0  0 0 1/2 0 0  0 0 1/2 0 0  0 0 -1 0 0];
  C1pz=C1p;
end
C2p=C1p'*C1p;
C2px=C1px'*C1px;
C2py=C1py'*C1py;
C2pz=C1pz'*C1pz;
C.Ptx=C1px;C.Pty=C1py;C.Ptz=C1pz;
%--- For projection--------
run cone;

g_bar=mean(reshape(g, nu*nv, na));

%-----------------------------------

%preparing for kronecker

Bx=construct_B(bspx,xi);
By=construct_B(bspy,yi);
Bz=construct_B(bspz,zi);
Bt=construct_B(bspt,ti);

[xxi,yyi,zzi]=ndgrid(xi,yi,zi);
%xxi=single(xxi);yyi=single(yyi);zzi=single(zzi);
vol=xsize*ysize*zsize;

%----single data type initialization------
geomx=single(zeros(size(ct)));
geomy=geomx; geomz=geomx;
img_out=geomx;
imggradx=geomx; imggrady=geomx; imggradz=geomx;

if(step==2)
    load temp1_new.mat;%temp884_lev2.mat;%results/dov16168_lev2_0601.mat;%temp16168_lev2;
    theta_x=2*xcoeff(:); theta_y=2*ycoeff(:); theta_z=2*zcoeff(:);
 %   load iter40_wo_ini.mat;
end;
