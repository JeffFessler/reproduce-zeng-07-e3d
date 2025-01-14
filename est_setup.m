%---------------------------------------------------------------
% Rigid setup difference estimation:[angx, angy, angz, tx, ty tz]
%
%Inputs
%   ctfile: the location of the CT volume
%   projfile: the location of the projection view file
%   angfile: the location of the projection angle file
%   idg: The indices of the projection views that correspond to deep-exhale
%       state.
%   ctfile_transformed: the location to store the setup difference
%       corrected CT volume.
%
%Output
%   img_out: transformed CT volumes, saved in the file specified by 
%       'ctfile_transformed'.
%User specified parameters
%   May need to tune the tradeoff parameters 'beta' and 'beta1'.
%   May change the threshhold 'ftol' for the termination condition.
%---------------------------------------------------------------
clear all;
ctfile='/n/ir24/y/rzeng/dov/phan/ct2.fld';
projfile='/n/ir24/y/rzeng/dov/phan/dynali_aperiod.fld';
angfile='/n/ir24/y/rzeng/dov/phan/dynaproj_aperiod_angs.fld';
ctfile_transformed= '/n/ir24/y/rzeng/dov/phan/ct_aff.fld';
idg=[20 40 60];
%load images and projection views
ct=fld_read(ctfile);
dsx=2;dsz=2;
ct=downsample3d(ct,dsx,dsz);
img_coeff=coeff_3d_mex(single(ct),3);
img_coeff=double(img_coeff);

%ct=downsample3d(ct,2,1);
proj=fld_read(projfile);%statprojli_ct2.fld');%%dynali_aperiod.fld');
projang=fld_read(angfile);%statproj_ct2_angs.fld');%
dsp=2;
g=proj(end:-1:1,:,idg);
g=downsample3d(g,dsp,1);
[ns,nt,np]=size(g);
wg=ones(ns,nt);
wg(:,round(1/2*nt):end)=0;
%wg(:,round(2/3*nt):end)=0;
maskg=logical(wg);

[nx,ny,nz]=size(ct);
xi=[0:nx-1];
yi=[0:ny-1];
zi=[0:nz-1];
[xxi,yyi,zzi]=ndgrid(xi,yi,zi);
img_grid=single([xxi(:) yyi(:) zzi(:)]');
clear xxi yyi zzi;

%generate Gb
%if(fp=='dd')  %Distance-driven forward projection method
%       orbit_start=182-(idg(1)-1)*358/666*4;
%       orbit=-(idg(end)-idg(1)+30)*358/666*4;
       % na=length(idg);
        orbit_start=projang(idg); 
        orbit_start=orbit_start(:);
        orbit=0;
      
        dx=0.9375*dsx; dy=0.9375*dsx; dz=3*dsz;
        ds=0.3880*2*dsp; dt=0.3880*2*dsp;
      
        [xsize,ysize,zsize]=size(ct);
        [ns,nt,na]=size(g);
        ig=image_geom('nx',xsize,'ny',ysize,'nz',zsize,'dx',dx,'dz', dz,...
            'offset_x',0,'offset_y',0,'offset_z',-0/dz);
        cg=ct_geom('fan','ns',ns,'nt',nt,'na',na,...
                    'ds',ds,'dt',dt,...
                    'dsd',1499.4,'dso',1000,...
                    'orbit', orbit,'orbit_start', orbit_start,...
                    'offset_s',-0.135/ds,'offset_t',-2.883/dt);
        G3=Gtomo_dd(cg,ig,'is_ns_nt',logical(1),'nthread',2);
   
    Gb = Gblock(G3,na);
    clear mask G3;
        
        
%         geo={'ns',ns,'nt',nt,'na',na,...
%         'dis_foc_src',  0,... % 0 for 3rd gen CT
%       'dis_src_det',1499.4,'dis_src_iso',1000,'dis_iso_det',499.4,...
%       'voxel_size',[dx dy dz],'ds', ds, 'dt', dt ,...
%       'channel_offset', -0.135/ds,... %s-direction detector offset
%       'offset_t',-2.883/dt,... %t-direction detector offset
%       'orbit', orbit, 'orbit_start', orbit_start,...
%      'img_offset',single([0 0 -8/dz]),...%-8/dz]),...
%       'is_ns_nt',logical(1)};
%     
%      mask=logical(ones([nx ny nz]));
%      G3=Gtomo3_dd(mask,geo{:});
%      Gb = Gblock(G3,na);
%    end;

%
iter=1;
angx=0; angy=0; angz=0; tx=0; ty=0; tz=0;
beta=2;
beta1=0.0001;
ftol=0.00002;
fdif=1; d_old=1;

while(abs(fdif)>ftol)
    tic
    %if(iter>50) break; end;
    rx=[1 0 0;0 cos(angx) -sin(angx); 0 sin(angx) cos(angx)];
    ry=[cos(angy) 0 sin(angy);0 1 0; -sin(angy) 0 cos(angy)];
    rz=[cos(angz) -sin(angz) 0; sin(angz) cos(angz) 0; 0 0 1];
    tr=[tx ty tz]';
    geo=rz*ry*rx*img_grid+kron(ones(1,nx*ny*nz),tr);
    
    r_grad_angx=rz*ry*[0 0 0; 0 -sin(angx) -cos(angx); 0 cos(angx) -sin(angx)];
    r_grad_angy=rx*[-sin(angy) 0 cos(angy); 0 0 0; -cos(angy) 0 -sin(angy)]*rz;
    r_grad_angz=rx*ry*[-sin(angz) -cos(angz) 0; cos(angz) -sin(angz) 0; 0 0 0];
    
    %image interpolation
    img_out=single(zeros(nx*ny*nz,1));
    imggradx=img_out; imggrady=img_out; imggradz=img_out;

    [img_out,imggradx, imggrady, imggradz]=...
        interpolate(img_coeff,double(geo(1,:)'),...
        double(geo(2,:)'),...
        double(geo(3,:)'), 3, 1,1,1);

    img_out=single(reshape(img_out,nx,ny,nz));
    img_out(find(img_out<0))=0;

    
    d=0; 
    grad_angx=0; grad_angy=0; grad_angz=0;
    grad_tx=0; grad_ty=0; grad_tz=0;
    %bp_gp = A'((p_m - p_bar)/T1 - (p_m - p_bar)/T2)
    for i=1:length(idg)
        g_mask=g(:,:,i).*wg;
        g_bar=mean(g_mask(maskg(:)));
        p=Gb{i}*img_out;
        p_mask=p.*wg;
        d=d - log(corr2(g_mask(maskg(:)),p_mask(maskg(:))));
        clear g_stack;
     
        p_bar = mean(p_mask(maskg(:)));
        T1 = sum(sum((g_mask(maskg(:))-g_bar).*(p_mask(maskg(:))-p_bar)));
        T2 = sum(sum((p_mask(maskg(:)) - p_bar).^2));  
        bp_gp = Gb{i}'*(((p_mask-p_bar).*wg/T2-(g_mask-g_bar).*wg/T1).*wg);
    %gradient
        img_grad_angx=imggrady'.*(r_grad_angx(2,:)*img_grid) + imggradz'.*(r_grad_angx(3,:)*img_grid);
        img_grad_angy=imggradx'.*(r_grad_angy(1,:)*img_grid) + imggradz'.*(r_grad_angy(3,:)*img_grid);
        img_grad_angz=imggradx'.*(r_grad_angz(1,:)*img_grid) + imggrady'.*(r_grad_angz(2,:)*img_grid);
        grad_angx=img_grad_angx*bp_gp(:) + grad_angx;
        grad_angy=img_grad_angy*bp_gp(:) + grad_angy;
        grad_angz=img_grad_angz*bp_gp(:) + grad_angz;
        grad_tx = imggradx'*bp_gp(:) + grad_tx;
        grad_ty = imggrady'*bp_gp(:) + grad_ty;
        grad_tz = imggradz'*bp_gp(:) + grad_tz;
        
        tmp=Gb{i}*imggradx;%img_grad_angx;
        hx = sum(tmp(maskg(:)).^2)/T2 - sum(tmp(maskg(:)))^2/sum(wg(:))/T2;
    end;

    angx = angx - beta1*grad_angx;
    angy = angy - beta1*grad_angy;
    angz = angz - beta1*grad_angz;
    tx = tx - beta*grad_tx;
    ty = ty - beta*grad_ty;
    tz = tz - beta*grad_tz;
    [angx angy angz tx ty tz]    
    dist(iter)=d;
    fdif = d - d_old
    d_old = d
    iter=iter+1
    toc
end;
 
fld_write(ctfile_transformed, img_out);
