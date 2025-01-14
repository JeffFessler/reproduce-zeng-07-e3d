function out=rm_air_tbl(in,thr,indx,indy)
% out=rm_air_tbl(in,thr,indx,indy)
%in: ct volume input
%thr: threshhold for rm_air
%indx: 2-element vector; leftmost and right most position of table
%indy: vector, y-axis search position for each slice

out = rm_air(in,thr);
zsize = size(in,3);
ysize = size(in,2);

%remove table
for k=1:zsize
    slice = out(:,:,k);
    for i=indx(1):indx(2)
        for j=indy(k):ysize
            %starting from indy(k), if the intensity value at some point 
            %is 0, then thorax region ends here so set the intensity value 
            %of the points below that point to be all 0. 
            if (slice(i,j)==0)
                slice(i,j:end)=0;
                break;
            end;
        end;
    end;
    out(:,:,k) = slice;
end;


