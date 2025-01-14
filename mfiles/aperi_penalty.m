function [pen, grad]=aperi_penalty(I,J,K,M,Pt, in,isgrad)

if(nargin==6) isgrad=1; end;
theta = reshape(in, I, J, K, M);
N = size(Pt,1);

pent = kron_product4(eye(I), eye(J), eye(K), Pt, in);

pen = sum(pent(:).^2)/(I*J*K*N);
grad=0;
if(isgrad==1)
    pen_gt = kron_product4(eye(I), eye(J), eye(K), Pt'*Pt, in)/(I*J*K*N)*2;
    grad = pen_gt(:);
end
