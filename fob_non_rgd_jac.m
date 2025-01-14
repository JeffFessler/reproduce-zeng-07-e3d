
%unmasked registration
function [f, img, geomx, geomy, geomz, pen,  grx, gry, grz, hess_x, hess_y, hess_z, jacobian]=fob_non_rgd_jac(theta_x, theta_y, theta_z, theta_x_locx, theta_x_locy, theta_x_locz, sub_image,coeff,lambda, reg_param, slop,sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, xi, yi, zi, jac_thr, lev)

%masked registration
%function [f, img, geomx, geomy, geomz, pen,  grx, gry, grz, hess_x, hess_y, hess_z, jacobian]=fob_non_rgd_jac(theta_x, theta_y, theta_z, theta_x_locx, theta_x_locy, theta_x_locz, sub_image,coeff,lambda, reg_param, slop,sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, xi, yi, zi, jac_thr, mask, range, lev)

%masked + affine
%function [f, img, geomx, geomy, geomz, pen,  grx, gry, grz, hess_x, hess_y, hess_z, jacobian, gax, gay, gaz, gt, hax, hay, haz, ht] = fob_non_rgd_jac( theta_x, theta_y, theta_z, theta_x_locx, theta_x_locy, theta_x_locz, sub_image,coeff,lambda, reg_param, slop,sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, xi, yi, zi, jac_thr, mask, range, lev ,ax, ay, az, t);

scale=2;


[x_size y_size z_size] = size(sub_image);
n_size = x_size*y_size*z_size;


%keyboard
%xo=x_size/2; yo=y_size/2; zo=z_size/2;

xcoeff=theta_x';
xloc_x=theta_x_locx;
xloc_y=theta_x_locy;
xloc_z=theta_x_locz;


ycoeff=theta_y';
yloc_x=theta_x_locx;
yloc_y=theta_x_locy;
yloc_z=theta_x_locz;

zcoeff=theta_z';
zloc_x=theta_x_locx;
zloc_y=theta_x_locy;
zloc_z=theta_x_locz;

%keyboard

sp=3;

max(max(max(theta_x)))

%For unmasked registration
%[mse img_trn img_gradx img_grady img_gradz x_deform y_deform z_deform obj_grad_x obj_grad_y obj_grad_z obj_hess_x obj_hess_y obj_hess_z penalty pen_grad_x pen_grad_y pen_grad_z pen_hess_x pen_hess_y pen_hess_z jacob]...
%    =grad_non_rgd_lim_jac_vol2(coeff, sub_image, xi,yi,zi,xcoeff,xloc_x,xloc_y,xloc_z,ycoeff,yloc_x,yloc_y,yloc_z,zcoeff,zloc_x,zloc_y,zloc_z,3, sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, reg_param, slop, jac_thr);%,mask,range);

%For multiresolution unmasked
[mse img_trn img_gradx img_grady img_gradz x_deform y_deform z_deform obj_grad_x obj_grad_y obj_grad_z obj_hess_x obj_hess_y obj_hess_z penalty pen_grad_x pen_grad_y pen_grad_z pen_hess_x pen_hess_y pen_hess_z jacob]...
    =grad_non_rgd_lim_jac(coeff, sub_image, xi,yi,zi,xcoeff,xloc_x,xloc_y,xloc_z,ycoeff,yloc_x,yloc_y,yloc_z,zcoeff,zloc_x,zloc_y,zloc_z,3, sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, reg_param, slop, jac_thr);%,mask,range);

%For maksed registration
%[mse img_trn img_gradx img_grady img_gradz x_deform y_deform z_deform obj_grad_x obj_grad_y obj_grad_z obj_hess_x obj_hess_y obj_hess_z penalty pen_grad_x pen_grad_y pen_grad_z pen_hess_x pen_hess_y pen_hess_z jacob]...
%    =grad_non_rgd_lim_jac_vol2_mask_rz(coeff, sub_image, xi,yi,zi,xcoeff,xloc_x,xloc_y,xloc_z,ycoeff,yloc_x,yloc_y,yloc_z,zcoeff,zloc_x,zloc_y,zloc_z,3, sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, reg_param, slop, jac_thr,mask,range,lev);

%For mask+affine
%[mse img_trn img_gradx img_grady img_gradz x_deform y_deform z_deform obj_grad_x obj_grad_y obj_grad_z obj_hess_x obj_hess_y obj_hess_z penalty pen_grad_x pen_grad_y pen_grad_z pen_hess_x pen_hess_y pen_hess_z jacob gax, gay, gaz, gt, hax, hay, haz, ht]...
%    =grad_non_rgd_lim_jac_vol2_aff_mask(coeff, sub_image, xi,yi,zi,xcoeff,xloc_x,xloc_y,xloc_z,ycoeff,yloc_x,yloc_y,yloc_z,zcoeff,zloc_x,zloc_y,zloc_z,3, sub_smp_xax,sub_smp_yax,sub_smp_zax, x_bnd, y_bnd, z_bnd, reg_param, slop, jac_thr,mask,range,lev, ax, ay, az, t);

mse =  mse

pen = penalty

f= mse + penalty

img = img_trn;


geomx = x_deform;
geomy = y_deform;
geomz = z_deform;

grx = obj_grad_x + pen_grad_x;
gry = obj_grad_y + pen_grad_y;
grz = obj_grad_z + pen_grad_z;
 
hess_x = obj_hess_x + pen_hess_x;
hess_y = obj_hess_y + pen_hess_y;
hess_z = obj_hess_z + pen_hess_z;

jacobian=jacob;
%keyboard
