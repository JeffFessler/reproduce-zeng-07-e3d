function f=cubic_bspline(x)
ax=abs(x);
if(ax<1)
    f=2/3-ax*ax+ax*ax*ax/2;
elseif(ax>=1)
    f=((2-ax)^3)/6;
else
    f=0;
end;
