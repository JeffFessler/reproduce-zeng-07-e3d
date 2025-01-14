function  hess = calc_apen_hessa(I,J,K,M,C, in)
%C.Ptx,y,z: for different aperiodicity penalties along x, y, z direction
hess=0;
parasize = I*J*K*M;
%N=size(C.Pt,1);
for i=1:3
    if(i==1) Pt=C.Ptx;
    elseif(i==2) Pt=C.Pty;
    else Pt=C.Ptz;
    end    
    N=size(Pt,1);
    d = reshape(in([1:parasize]+(i-1)*parasize), I, J, K, M);
    ht = kron_product4(eye(I), eye(J), eye(K), Pt'*Pt, d)/(I*J*K*N)*2;
    hess = hess+d(:)'*(ht(:));
end;  
  

