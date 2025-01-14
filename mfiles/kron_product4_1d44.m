function out=kron_product4_1d44(b1,b2,b3,b4,A)
%out=kron_product4_1d44(b1,b2,b3,b4,A)
%Implement (b4@b3@b2@b1)*A, where @ stands for kronecker operator
%b4 is a 1D vector
%It is faster than kron_product4 when b4 contains many zeros
%and overcome the possible memory blow-up problem
%Sep. 2, by rzeng

%clear d;
if (size(b4,1)>1 & size(b4,2)>1)
    disp 'The 4th parameter must be 1D';
    return;
end;
n1=size(b1,1); n2=size(b2,1); n3=size(b3,1); n4=size(b4,1);
l1=size(b1,2); l2=size(b2,2); l3=size(b3,2); l4=size(b4,2);
B=reshape(A,l1,l2,l3,l4);
ind=find(b4);
b4new=b4(ind);
for j=1:length(ind)
    for i=1:l3
        b=B(:,:,i,ind(j));
        c1(:,i)=reshape(b1*b*b2',n1*n2,1);
    end
    d(:,j)=reshape(c1*b3',n1*n2*n3,1);
end

if (isempty(ind))
    out=single(zeros(n1,n2,n3,n4));
else
    c2=d*b4new'; 
    out=reshape(c2,n1,n2,n3,n4);        
end;