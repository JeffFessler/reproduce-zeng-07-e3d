%sig
function cyc=find_breathcyc(est)
if(~exist('est'))
    est=resp_proj;
end
h=[1 1 1 1 1 1 1]/7;
estf=conv(h,est);
estf=estf(4:end-3);
%estf(1)=estf(2);
%estf(end)=estf(end-1);


sig=estf;
peak=[ ];
valley=[ ];
for i=2:(length(sig)-1)
    if ( sig(i)>=sig(i-1) & sig(i)>=sig(i+1) )
        peak=[peak i];
    end;
    if( sig(i)<=sig(i-1) & sig(i)<=sig(i+1) )
        valley=[valley i];
    end;
end;

figure;
plot(estf);title 'Breathing signal')
disp 'hit any key to continue');
pause
r=input('how many breathing cycles in the breathing signal:\n','s')
ncyc=str2num(r);
%remove invalid peaks
%ncyc=4;
[temp,J]=sort(sig(peak),'descend');
peak=sort(peak(J(1:ncyc)));
clear temp
%[temp,J]=sort(sig(valley),'ascend');
%valle=sort(valley(J(1:ncyc)));
%valley=valle;
%id=[ ];
%for i=1:length(peak)
%    if sig(peak(i))>=avgpeak
%        id=[id i];
%    end;
%end;
%peak=peak(id);

%between each two peaks, pick one valley
j=1;
for i=1:length(peak)-1
    idx=find(valley>peak(i) & valley<peak(i+1));
    n=length(idx);
    if n>1
        val(i)=find(sig==min(sig(valley(j:j+n-1))));
        j=j+n;
    else
        val(i)=valley(j);
        j=j+1;
    end;
end;

scale=20; % below 5% represents exhale
ncyc=length(peak);
cyc=zeros(ncyc,2);
cyc(1,1)=0; cyc(end,end)=length(estf)-1;
for i=1:ncyc-1
    x=sig(peak(i):peak(i+1));
    ex1=(sig(peak(i))-sig(val(i)))/scale;
    ex2=(sig(peak(i+1))-sig(val(i)))/scale;
    idx1=find(x<sig(val(i))+ex1);
    idx2=find(x<sig(val(i))+ex2);
    id1=find(x==sig(val(i)));
    id2=find(x==sig(val(i)));
    cyc(i,2)=val(i)-(id1-idx1(1));
    cyc(i+1,1)=val(i)+(idx2(end)-id2);
end;
x=sig(1:peak(1));
ex=(sig(peak(1)) - sig(1))/scale;
idx=find(x<sig(1)+ex);
cyc(1,1)=idx(end);

x=sig(peak(end):end);
ex=(sig(peak(end)) - sig(end))/scale;
idx=find(x<sig(end)+ex);
cyc(end,2)=length(sig)-(length(x)-idx(1));

cyc=cyc-1; %t index starts from 0 rather than 1 in matlab array