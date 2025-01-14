%-------------------------------------------------
%Convert the projection views from photon counts to attenuation.
%
%Inputs:
%   prjfile: the location of the raw projection file
%   tblfile: the location of the raw table scan projection file
%   attfile: the location to save the DRRs.
%
%-----------------------------------------------------
prjfile='/n/ir24/y/rzeng/dov/temp/phan/dynaproj_aperiod.fld';
tblfile='/n/ir24/y/rzeng/dov/temp/phan/proj_tbl.fld';
drrfile='/n/ir23/y/rzeng/dov/temp/phan/dynali_aperiod.fld';
proj=fld_read(tblfile);
tbl=single(proj(5:end-4,5:end-4,1:2:end));
proj=fld_read(prjfile);
prj=single(proj(5:end-4,5:end-4,:));

prjwotbl = log(tbl ./ prj);
fld_write(drrfile, prjwotbl);

