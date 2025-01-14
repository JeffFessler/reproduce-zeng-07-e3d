function  hess = calc_rpen_hessa(Px, Py, Pz, Pt, in)
hess=0;
I=size(Px,2);
J=size(Py,2);
K=size(Pz,2);
M=size(Pt,2);
I1=size(Px,1);
J1=size(Py,1);
K1=size(Pz,1);
M1=size(Pt,1);

parasize = I*J*K*M;
for i=1:3
    d = reshape(in([1:parasize]+(i-1)*parasize), I, J, K, M);
    hx = kron_product4(Px'*Px, eye(J), eye(K), eye(M), d)/(I1*J*K*M)*2;
    hy = kron_product4(eye(I), Py'*Py, eye(K), eye(M), d)/(I*J1*K*M)*2;
    hz = kron_product4(eye(I), eye(J), Pz'*Pz, eye(M), d)/(I*J*K1*M)*2;
    ht = kron_product4(eye(I), eye(J), eye(K), Pt'*Pt, d)/(I*J*K*M1)*2;
    hess = hess+d(:)'*(hx(:)+hy(:)+hz(:)+ht(:));
end;  
  

