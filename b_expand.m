% length n -> 2n-1
  %  load temp3_60.mat;%result_lev2_new;
    xl=length(y_locx);
    yl=length(y_locy);
    zl=length(y_locz);

    theta_x_init = reshape(theta_x,yl,xl,zl);
    theta_x_init = transpose_3D(theta_x_init);
    theta_y_init = reshape(theta_y,yl,xl,zl);
    theta_y_init = transpose_3D(theta_y_init);
    theta_z_init = reshape(theta_z,yl,xl,zl);
    theta_z_init = transpose_3D(theta_z_init);
    
    new_thetax=2*upsample_by_2(theta_x_init,3);
    new_thetay=2*upsample_by_2(theta_y_init,3);
    new_thetaz=2*upsample_by_2(theta_z_init,3);
    
    new_thetax = transpose_3D(new_thetax);
    xcoeff=new_thetax(:);
    new_thetay = transpose_3D(new_thetay);
    ycoeff=new_thetay(:);
    new_thetaz = transpose_3D(new_thetaz);
    zcoeff=new_thetaz(:);
    
    new_locx=zeros(1,2*xl-1);
    new_locy=zeros(1,2*yl-1);
    new_locz=zeros(1,2*zl-1);
    new_locx(1:2:2*xl-1)=y_locx*2; new_locx(2:2:2*xl-2)=[y_locx(1:end-1)+y_locx(2:end)];
    new_locy(1:2:2*yl-1)=y_locy*2; new_locy(2:2:2*yl-2)=[y_locy(1:end-1)+y_locy(2:end)];
    new_locz(1:2:2*zl-1)=y_locz*2; new_locz(2:2:2*zl-2)=[y_locz(1:end-1)+y_locz(2:end)];
    
    if(step==1)
    eval(['save ' tmpfile '/temp1_new xcoeff ycoeff zcoeff new_locx new_locy new_locz']);
    end;
    if(step==2)
    eval(['save ' tmpfile '/temp2_new xcoeff ycoeff zcoeff new_locx new_locy new_locz']);
    end;
    %save bsp_coeff_lev2_new xcoeff ycoeff zcoeff;
    disp "Finished!";
