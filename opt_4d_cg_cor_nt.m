%-----------------------
%DOV with B-spline motion model: Correlation-base metric, CG algorithm with
%one step newton method to find the stepsize instead of linesearch
%Date: 05/17/2006
%E=sum_over_m(-(log(cor(g_m,p_m)))+penalty
%   g_m: measured projection views
%   p_m: modelled projection views
%
%Parameters that may be tuned by users
%   tolf, alpha, beta1, beta
%-----------------------

%clear all;

dovinit=1; step=2;
run setup_w_ini;

%---optimization--------------------
tolf=1e-5;%1e-4;%1e-5;%1e-6;%10;

alpha = 4;%0.5;%3*1e-6;%0.05;%%500;

beta1=0.0001; beta=0.00001;
if (dovinit==1) beta=2e-5; beta1=0.001;  end;

thr = 7;
fold=100;
fdif=tolf+1;
dif_theta=1;
iter=0;
pen=0;
iter=1;
d1=3000;

count=0; %count for update times
update=0;
maxxgrad=1;
ob_bsp.Bx=Bx;ob_bsp.By=By;ob_bsp.Bz=Bz;ob_bsp.Bt=Bt;
ndiv=0; %for inexact line search


%n=round(ndiv/2);
while (abs(fdif)>tolf) %(dif_theta>0.03) %(d1>2100)%
%    if iter==2 load theta_noapp; end
     xgrad = zeros(parasize,1);
     ygrad = xgrad;
     zgrad = xgrad;
     
    iter


 %  profile on

d1=0; 
tic
    %caiculate the CC and the derivatives
    [d1,xgrad,ygrad,zgrad]=calc_corr_grad(ob_bsp, theta_x, theta_y, theta_z, g, g_bar,...
                    img_coeff,xxi,yyi,zzi, Gb, 1);


    
    disp('The cost is :..............');
    dd(iter)=d1;
   
    %------Calculate the penalties and their derivatives -----------------------------
    [rpenx, rpen_xgrad] = rough_penalty(C1x, C1y, C1z, C1t, theta_x);
    [rpeny, rpen_ygrad] = rough_penalty(C1x, C1y, C1z, C1t, theta_y);
    [rpenz, rpen_zgrad] = rough_penalty(C1x, C1y, C1z, C1t, theta_z);
    
    [apenx, apen_xgrad] = aperi_penalty(lx,ly,lz,ltt, C1px, theta_x);
    [apeny, apen_ygrad] = aperi_penalty(lx,ly,lz,ltt, C1py, theta_y);
    [apenz, apen_zgrad] = aperi_penalty(lx,ly,lz,ltt, C1pz, theta_z);
    pen1 = beta*(rpenx +rpeny +rpenz)
    pen2 = beta1*(apenx+apeny+apenz)
    pen = pen1+pen2; %beta*(rpenx+rpeny+rpenz)+ beta1*(apenx+apeny+apenz)
    pen_xgrad = beta*rpen_xgrad + beta1*apen_xgrad; 
    pen_ygrad = beta*rpen_ygrad + beta1*apen_ygrad; 
    pen_zgrad = beta*rpen_zgrad + beta1*apen_zgrad;

    dist(iter)=d1+pen;
    dist(iter)
    fdif(iter)=dist(iter)-fold;
    %fold=dist(iter); 
toc
 
if(fdif(iter)>0) disp 'reset'; reset=1; fold=dist(iter); continue; end;

 tic   
    theta = [theta_x;theta_y;theta_z];
     grad = [xgrad(:); ygrad(:); zgrad(:)]+ [pen_xgrad(:); pen_ygrad(:); pen_zgrad(:)];
     
    %Conjugate gradient algorithm
     if (iter==1)
         r=0; d_old=zeros(3*parasize,1);
     else
         y = grad-grad_old;
         r=real(grad'*(y-2*d_old*(y'*y)/(d'*y))) / real(d_old'*y);
     end
     d = -grad + r*d_old;             
         
     theta_old=theta;
     
     %stepsize
     grad_a = d'*grad;
     L_hess_a = calc_cor_hess_a(ob_bsp, theta_x, theta_y, theta_z,...
                  img_coeff,na, xxi,yyi,zzi, Gb, d);
     P_hess_a = beta* calc_rpen_hessa(C1x, C1y, C1z, C1t, d) + ...
           beta1* calc_apen_hessa(lx,ly,lz,ltt,C, d);
     alpha0 = -(grad_a)/(L_hess_a+P_hess_a)
        
     dir = alpha0*d;     
     theta = theta_old + dir;
     
        theta_x = theta(1:parasize);
        theta_y = theta([1:parasize]+parasize);
        theta_z = theta([1:parasize]+2*parasize);       
        grad_old=grad;
        d_old = d;
        fold=dist(iter)
        dif_theta=max(abs(theta-theta_old))
    toc    
  %   profile viewer
  
    save temp theta_x theta_y theta_z iter;
     iter=iter+1;
     beep
end;
save temp.mat theta_x theta_y theta_z;
%save /n/ir24/y/rzeng/temp.mat theta_x theta_y theta_z img_out geomx geomy geomz;
beta
beta1
beep;
%save spl_xyzt.mat dist d1 pen corr theta iter geomx geomy xcoeff_ini ycoeff_ini;


clear partial_dirx partial_diry partial_dirz partial_geomx partial_geomy partial_geomz;

if(step==1)
  xcoeff=bsp_expand(theta_x,[lx ly lz ltt],[1 1 1 0]);
  ycoeff=bsp_expand(theta_y,[lx ly lz ltt],[1 1 1 0]);
  zcoeff=bsp_expand(theta_z,[lx ly lz ltt],[1 1 1 0]);
  save temp1_new.mat xcoeff ycoeff zcoeff;
  disp 'temp1_new.mat saved'
end

if(step==2)
    save test theta_x theta_y theta_z img_coeff xi yi zi Bx By Bz Bt Gb;
end;
%complete=step;
