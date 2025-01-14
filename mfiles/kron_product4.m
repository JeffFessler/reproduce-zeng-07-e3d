function out=kron_product4(b1,b2,b3,b4,A)
%out=kron_product4(b1,b2,b3,b4,A)
%Implement (b4@b3@b2@b1)*A, where @ stands for kronecker operator
%Apr.05, by rzeng
clear d;
n1=size(b1,1); n2=size(b2,1); n3=size(b3,1); n4=size(b4,1);
l1=size(b1,2); l2=size(b2,2); l3=size(b3,2); l4=size(b4,2);
B=reshape(A,l1,l2,l3,l4);
for j=1:l4
    for i=1:l3
        b=B(:,:,i,j);
        c1(:,i)=reshape(b1*b*b2',n1*n2,1);
    end
    d(:,j)=reshape(c1*b3',n1*n2*n3,1);
end
c2=d*b4';
 
out=reshape(c2,n1,n2,n3,n4);        