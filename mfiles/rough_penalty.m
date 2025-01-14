function [pen, grad]=rough_penalty(Px, Py, Pz, Pt, in, isgrad)

if(nargin==5) isgrad=1; end;

I=size(Px,2);
J=size(Py,2);
K=size(Pz,2);
M=size(Pt,2);
I1=size(Px,1);
J1=size(Py,1);
K1=size(Pz,1);
M1=size(Pt,1);

theta = reshape(in, I, J, K, M);
penx = kron_product4(Px, eye(J), eye(K), eye(M), in);
peny = kron_product4(eye(I), Py, eye(K), eye(M), in);
penz = kron_product4(eye(I), eye(J), Pz, eye(M), in);
pent = kron_product4(eye(I), eye(J), eye(K), Pt, in);

pen = sum(penx(:).^2/(I1*J*K*M) + peny(:).^2/(I*J1*K*M) + penz(:).^2/(I*J*K1*M) + pent(:).^2/(I*J*K*M1));

grad=0;
if(isgrad==1)
pen_gx = kron_product4(Px'*Px, eye(J), eye(K), eye(M), in)/(I1*J*K*M)*2;
pen_gy = kron_product4(eye(I), Py'*Py, eye(K), eye(M), in)/(I*J1*K*M)*2;
pen_gz = kron_product4(eye(I), eye(J), Pz'*Pz, eye(M), in)/(I*J*K1*M)*2;
pen_gt = kron_product4(eye(I), eye(J), eye(K), Pt'*Pt, in)/(I*J*K*M1)*2;
grad = pen_gx(:) + pen_gy(:) + pen_gz(:) + pen_gt(:);
end;
