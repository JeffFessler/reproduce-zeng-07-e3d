function out = calc_fp(A,f,i)
%To calculate: out = A{i}f
%A: system matrix
%f: deformed ct volumes
%i: projection number
[nx,ny,nz] = size(f);
vol=nx*ny*nz;

%out = A{i}*reshape(f,vol,1); %3l
out = A{i}*f; %dd

out = single(out);
