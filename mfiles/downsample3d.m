function out=downsample3d(in,dsx,dsz)
%Downsample by averaging
%input: 
%    in: input array
%    dsx: tranxasial downsampling by averaging
%    dsz: z direction decimating by a factor of 'dsz'
%
%By rzeng, 2005
[xsize,ysize,zsize]=size(in);
for i=1:dsz:zsize
    j=(i-1)/dsz+1;
    out(:,:,j)=downsample2(in(:,:,i),dsx);
end;
