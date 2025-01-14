function B = construct_B(knot, xi)
%function B = constru_B(knot, xi)
%Construct a Bspline matrix.
%knot is a structure variable: 
%             knot.loc: a vector containing the locations of Bspline knots 
%             knot.h: a vector defining the width of each Bspline function
%xi: a vector giving the resampling locations of the signal represented by
%    the Bsplines.

nknot=length(knot.loc);
nx=length(xi);
B=zeros(nx,nknot);

for i=1:nx
    for j=1:nknot
        x=(xi(i)-knot.loc(j))/knot.h(j);
        if(abs(x)<2)
            B(i,j)=cubic_bspline(x);
        end
    end
end
