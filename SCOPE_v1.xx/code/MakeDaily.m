function [d,quality,days] = MakeDaily(time,data,m,QC,Nmin)

%function d=MakeDaily(time,data,start,finish,QC,Nmin) converts a matrix of
%data into daily average (m=0) or sum values (m=1), minimum (m=2) or
%   maximum (m=3)
%The output: daynumber
%next columns: daily data
%last column:zeros and ones, where 1 means that there were missing
%values during this day.
%optional: give start and finish (integer julian daynumbers)
%optional: give QC (quality control): two numbers indicating the minimum and maximum value that can be accepted
% Nmin: the minimum number of data points on a day

start = ceil(min(time));
finish = floor(max(time));

%interval = time(2)-time(1);
s2 = size(data,2);
if nargin<5
    Nmin = 1;
end

[d,quality] = deal(NaN*zeros(finish-start+1,s2));
days = (start:finish)';
for i = 1: finish-start+1
    for j=1:s2
        if nargin>3
            f=find((time<=i+start) & (time>i-1+start) & data(:,j)<QC(2) & data(:,j)>QC(1));
        else
            f=find((time<=i+start) & (time>i-1+start));
        end
        z = find(~isnan(data(f,j)));
        if length(z)>Nmin
            
            switch m
                case 0,
                    d(i,j)=mean(data(f(z),j));
                case 1,
                    d(i,j)=sum(data(f(z),j));
                case 2,
                    d(i,j)=min(data(f(z),j));
                otherwise,
                    d(i,j)=max(data(f(z),j));
            end
        else
            d(i,j) = NaN;
        end   
    end
    quality(i) = length(z);
end
