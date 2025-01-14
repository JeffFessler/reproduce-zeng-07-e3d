%----------------
%LS fitting of the simplified motion estimate into B-Spline motion
%
%Inputs
%   geomx, geomy, geomz: extreme deformation fields from image
%       registration.
%   a: the proportionality sequence estimated from 'dovinit_cor.m'
%
%Output
%   theta_x_ini, theta_y_ini, theta_z_ini: the bspline coeffcients
%   bspx, bspy, bspz, bspt: the Bspline knot locations
%   cyc, knotpercyc: the number of breathing cycles in the scanning period
%       the number of temporal knots in one breathing cycle
%   (the outputs are saved in 'dovinit.mat'.)
%
%User specified parameters
%   hx, hy, hz: knot spacings along x, y, z directions respectively
%----------------
 if(1)
 load data/reg_ct02_ds221_16168.mat; %deformation fields from registration
  geomx=single(geomx);
  geomy=single(geomy);
  geomz=single(geomz);
  load data/a.mat;     % estimated proportionality sequence
  
  %cyc=[0 19; 22 35; 39 58; 63 77];
  cyc=find_breathcyc(a);%[0 19; 23 36; 40 60; 63 79];
  knotpercyc=5;
  
  [nx,ny,nz]=size(geomx);
  nt=length(a);
  
  %Bspline knot placement
  hx=16; hy=16; hz=8;
  xi=linspace(0,nx-1,nx); yi=linspace(0,ny-1,ny); zi=linspace(0,nz-1,nz);ti=linspace(0,nt-1,nt); 
  
  bspx.loc=[0:hx:nx-1]; bspy.loc=[0:hy:ny-1]; bspz.loc=[nz-1:-hz:0]; 
  bspx.h=hx*ones(size(bspx.loc));bspy.h=hy*ones(size(bspy.loc));bspz.h=hz*ones(size(bspz.loc));
 
  nt=length(a);
  bspt.loc = [linspace(cyc(1,1),cyc(1,2),knotpercyc) linspace(cyc(2,1),cyc(2,2),knotpercyc) ...
          linspace(cyc(3,1),cyc(3,2),knotpercyc) linspace(cyc(4,1),cyc(4,2),knotpercyc)];
  bspt.h = [ones(1,knotpercyc)*(cyc(1,2)-cyc(1,1))/(knotpercyc-1) ones(1,knotpercyc)*(cyc(2,2)-cyc(2,1))/(knotpercyc-1) ...
         ones(1,knotpercyc)*(cyc(3,2)-cyc(3,1))/(knotpercyc-1) ones(1,knotpercyc)*(cyc(4,2)-cyc(4,1))/(knotpercyc-1)]; 
  
 bsp_xyz{1}=bspx;
 bsp_xyz{2}=bspy;
 bsp_xyz{3}=bspz;
 bsp_t{1}=bspt;  

 %Fitting
 cx=bsp_fit(geomx,bsp_xyz);
 cy=bsp_fit(geomy,bsp_xyz);
 cz=bsp_fit(geomz,bsp_xyz);
 ct=bsp_fit(a,bsp_t);
 
 theta_x_ini=kron(ct,cx);
 theta_y_ini=kron(ct,cy);
 theta_z_ini=kron(ct,cz); 
  
 save ~rzeng/dov/data/dovinit.mat theta_x_ini theta_y_ini theta_z_ini bspx bspy bspz bspt...
      cyc knotpercyc;
 end