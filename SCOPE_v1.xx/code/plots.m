% plots.m (script) makes plots the output of SCOPE_v1.51 of the latest run.

clc, clear all, close all

directories         =   dir('..\output\*');
[time_value_s,I]    =   sort([directories(3:end).datenum]);
Directory           =   directories(2+I(end)).name;

%% load verification data
path1_ = ['..\output\' Directory ,'\'];
info1   = dir([path1_ '\*.dat']);           %the most recent output

L = length(info1);
wl = dlmread([path1_ 'wl.dat'],'',2,0);

for i = 1:L-1
    s1 = info1(1).bytes;
    n1 = info1(1).name;
    if ~ (strcmp(info1(i).name,'pars_and_input.dat') || strcmp(info1(i).name,'pars_and_input_short.dat'))
        
        D1 = dlmread([path1_ info1(i).name],'',2,0);
        
        spn = ceil(sqrt(size(D1,2)));
        h1 = textread([path1_ info1(i).name],'%s');
        figure(i)
        if spn>7
            nl = length(D1);
            for z = 1:min(47,nl)
                if size(D1,1)-1==length(wl)
                    plot(wl,D1(z,1:end-1)'), hold on
                    set(gca,'xlim',[.4 2.5])
                else
                    if z<size(D1,1)+1
                        plot(D1(z,:)','r'), hold on
                    end
                end
            end
            title(info1(i).name,'Interpreter', 'none')
        else
            for m = 1:size(D1,2)
                subplot(spn,spn,m)
                plot(D1(:,m),'r')
                title([info1(i).name h1(m)],'Interpreter', 'none')
            end
        end
        differentcontent = 1;
    end
end