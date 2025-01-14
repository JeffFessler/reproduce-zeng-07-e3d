fp = 'newdd';
    if(streq(fp,'3l'))  %3l forward projector, siddon's alg
        pixel_size=0.1956*10;
        nu=180; nv=200;
        su=2;sv=2; sz=0.2517/0.0978;
        orbit=180; orbit_start=0;
        systype=sprintf('3l@1000,500,%d,%d,%2.1f,%2.1f,%2.1f,0,0,0,0,0@-@-2d,%d,%d,%d', nu, nv,su,sv,sz,na,orbit, orbit_start);
        prj_mask = logical(ones([xsize ysize zsize]));
        G3=Gtomo3(systype, prj_mask, xsize, ysize, zsize);
        G3.cascade_after = pixel_size;
        tt = sprintf('%g * permute(x, [2 1 3])', G3.cascade_after);
	    G3.cascade_after = inline(tt, 'x', 'is_transpose', 'istart', 'nblock');
    end;
    if(streq(fp,'dd'))  %Distance-driven forward projection method
 %        dsx=2;dsz=2;dsp=1;
 %     nu=504; nv=376; na=84;
       [nu,nv,na]=size(g);
        dx=0.9375*dsx; dy=0.9375*dsx; dz=3*dsz;
        ds=0.3880*2*dsp; dt=0.3880*2*dsp;
        orbit_start=182; orbit=-na*4*358/666;
       % offx==1.7685; offy=-0.4999; offz=-7.0125;
        offx=0; offy=0; offz=0;
        geo={'ns',nu,...
        'nt',nv,...   %row number
        'na',na,...
        'dis_foc_src',  0,... % 0 for 3rd gen CT
      'dis_src_det',1499.4,...
      'dis_src_iso',1000,...
      'dis_iso_det',499.4,... % = arg.dis_src_det - arg.dis_src_iso;
      'voxel_size',[dx dy dz],...
      'ds', ds,...
      'dt', dt ,...%= 0.625 * arg.dis_src_det / arg.dis_src_iso;
      'channel_offset', -0.135/ds,... %s-direction detector offset
      'offset_t',-2.883/dt,... %t-direction detector offset
      'orbit', orbit,...     % + clockwise; - counter clockwise
      'orbit_start', orbit_start,...
     'img_offset',single([offx/dx offy/dy offz/dz]),... %single([0 0 -8/dz]),...
        'is_ns_nt',logical(1)};
    
     mask=logical(ones([xsize ysize zsize]));
     G3=Gtomo3_dd(mask,geo{:});
     Gb = Gblock(G3,na);
    end;
    
    if(streq(fp,'newdd'))%       
  %      dsx=2;dsz=2;dsp=1;
  %    nu=504; nv=376; na=167;
        
        dx=0.9375*dsx; dy=0.9375*dsx; dz=3*dsz; 
        ds=0.3880*2*dsp; dt=0.3880*2*dsp;
        [xsize,ysize,zsize]=size(ct);
        [ns,nt,na]=size(g);
        orbit=0;
        orbit_start=fld_read('/n/ir24/y/rzeng/dov/phan/dynaproj_aperiod_angs.fld');
        orbit_start=orbit_start(:);
        ig=image_geom('nx',xsize,'ny',ysize,'nz',zsize,'dx',dx,'dz', dz,...
            'offset_x',0,'offset_y',0,'offset_z',0);
        cg=ct_geom('fan','ns',nu,'nt',nt,'na',na,...
                    'ds',ds,'dt',dt,...
                    'dsd',1499.4,'dso',1000,...
                    'orbit', orbit,'orbit_start', orbit_start,...
                    'offset_s',-0.135/ds,'offset_t',-2.883/dt);
        G3=Gtomo_dd(cg,ig,'is_ns_nt',logical(1),'nthread',2);
    end;   
    
    Gb = Gblock(G3,na);
    clear mask G3;
    
%Check the CBCT geometry    
%     rgb=zeros(504,376,3);
% rgb(:,:,1)=p40/max(p40(:));
% rgb(:,:,2)=y(:,:,40)/max(y(:));
