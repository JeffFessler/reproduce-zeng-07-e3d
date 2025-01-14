%----------------------------------------------------
%Read the projection views and their projection angles into two files
%respectively.
%
%Inputs
%   str: the directory contains all the projectio views
%   prjfile: the location to store the projection view matrix
%   anggile: the location to store the projection angle vector
%
%Use sepcified parameters
%   ds: temporal downsampling rate of projection views
%-------------------------------------------------------
clear proj dir;
str='/n/ir24/y/rzeng/dov/newDOVexperiment/projections/2006.03.17.20_09_12/';
                                  %static proj. of ct2:2006.03.17.19_38_34
                                  %dynamic proj. of aperiodic motion:2006.03.17.20_09_12
prjfile='/n/ir24/y/rzeng/dov/phan/dynaproj_aperiod.fld';
angfile='/n/ir24/y/rzeng/dov/temp/phan/dynaproj_aperiod_angs.fld';
D = dir(str);
M=length(D);
disp('\n Begin, please wait...\n');
ds=4;%4;
n=floor((M-3)/2); 
        %(M-3) %for table;%(M-4)/2;
for i=1:ds:n
    tic
 %   disp(i);
    file=[str D(i+2).name];
 %read the image   
    p=hnd_read(file);
    proj(:,:,(i-1)/ds+1)=downsample2(single(p.pixels),2);
 
 %read the header
    p=hnd_readinfo(file);
    ang((i-1)/ds+1)=p.MachineCTProjectionAngle;
    clear p file;
    toc
end;
fld_write(prjfile,proj);%dynaproj_aperiod.fld', proj);

%ang=ang+90;
fld_write(angfile,ang);