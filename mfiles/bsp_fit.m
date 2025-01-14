function [theta,B] = bsp_fit(sig, bspgrid);
%function theta = bsp_fit(sig, bspgrid);
%bspgrid: cell structure
%   bspgrid{1}.loc, bspgrid{1}.h: x-axis direction spline knot distribution
%   bspgrid{2}.loc, bspgrid{2}.h: y-axis direction spline knot distribution
%   ...
%06/05/2006, by rzeng

dim=size(sig);
ndim=ndims(sig);
if(dim(1)==1|dim(2)==1) ndim=1; end;
if(ndim>4) 
    disp 'Not implemented fitting for signals with dimension>4!';
    theta=0;
    return;
end;
for i=1:ndim
    sig_grid=linspace(0,dim(i)-1,dim(i));
    B{i}=construct_B(bspgrid{i},sig_grid);
    BTBinv{i}=inv(B{i}'*B{i});
end;

switch(ndim)
    case {1}
        c=B{i}'*sig(:);%kron_product(Bx',By',Bz',Bt',gt);
        theta = BTBinv{1}*c;
    case {2}
        c=kron(B{1}',B{2}')*sig(:);
        theta = kron(BTBinv{1},BTBinv{2})*c;
    case {3}
        c=kron_product4(B{1}',B{2}',B{3}',1,sig);
        theta = kron_product4(BTBinv{1},BTBinv{2},BTBinv{3},1,c);
    case {4}
        c=kron_product4(B{1}',B{2}',B{3}',B{4}',sig);
        theta = kron_product4(BTBinv{1},BTBinv{2},BTBinv{3},BTBinv{4},c);
    otherwise
        disp 'Not implemented!'
        theta=0;
end;
