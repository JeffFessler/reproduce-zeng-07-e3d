function out = HU2mu(in,mu_water)
%out = HU2mu(in,mu_water)
if (max(in(:))<10)
    disp 'already in mu!'
end;
%mu_water = 0.2;
out = (in-1000)*mu_water/1000 + mu_water;