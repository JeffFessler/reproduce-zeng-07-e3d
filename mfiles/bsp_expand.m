function out = bsp_expand(theta,sz,dir)  
% out = bsp_expand(in,size,dir): expand bspline coefficients by 2.
% size=[xl,yl,zl,tl]: Size of the control knots
% dir=[dirx, diry,dirz,dirt]: dirx/y/z/t=0 or 1
%                             0 means no expansion
%                             1 means expansion
%Feb. 07, 2005
    xl=sz(1);
    yl=sz(2);
    zl=sz(3);
    tl=sz(4);
    
    dirx=dir(1);
    diry=dir(2);
    dirz=dir(3);
    dirt=dir(4);

    theta_new=[];
if(dirx==1 & diry==1 & dirz==1 & dirt==0)
    theta_init = reshape(theta,yl,xl,zl,tl);
    for i=1:tl
       theta_t = theta_init(:,:,:,i);
       theta_t = transpose_3D(theta_t);
       theta_t_new = upsample_by_2(theta_t,3);
       tmp = transpose_3D(theta_t_new);
       theta_new=[theta_new; tmp(:)];
   end
elseif(dirx==1 & diry==1 & dirz==0 & dirt==0)
    theta_init = reshape(theta,xl,yl,zl,tl);
    for i=1:tl
	for j=1:zl
       	    theta_zt = theta_init(:,:,j,i)';
            theta_zt_new = upsample_by_2(theta_zt,2);
            tmp = theta_zt_new';
            theta_new(:,:,j,i) = tmp;
	end
   end
else
       disp 'Not implement yet!'
end
   out=theta_new;
