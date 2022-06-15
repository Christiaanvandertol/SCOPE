function output_verification_csv(Output_dir, verification_dir)
% Date: 07 August 2012
% Author: Christiaan van der Tol (c.vandertol@utwente.nl)
% output_verification.m (script) checks if the output of the latest run
% with SCOPE_v1.51 matches with a 'standard' output located in a directory
% called 'verificationdata'. If it does not, warnings will appear in the
% Matlab command window.
% The following is tested:
%   - does the number of output files match?
%   - does the size of the files match (number of bytes)?
%   - are all files that are in the verification dataset present with the
%   same file names?
%   - is the content of the files exactly the same?
% If the output is different, for example because different parameter
% values have been used in the simulations, then the variables that are
% different will be plotted: the verification data in blue, and the latest
% run in red. In this way the differences can be visually inspected.

% clc, close all
%
% directories         =   dir(['..' filesep 'output' filesep '*']);
% [time_value_s,I]    =   sort([directories(3:end).datenum]);
% Directory           =   directories(2+I(end)).name;
%
% Directory = Output_dir

%% load verification data
path0_ = ['output' filesep verification_dir filesep];
path1_ = [Output_dir filesep];

info0   = dir([path0_ filesep '*.csv']);         %'standard' validation data (to compare with)
info1   = dir([path1_ filesep '*.csv']);           %the most recent output

[differentsize,differentcontent,differentnumberoffiles]  = deal(0);

if ~(length(info0)==length(info1))
    fprintf(['\nWarning: in the output file, ' num2str(length(info1)) ' files were stored, \r'])
    fprintf(['whereas there should be ' num2str(length(info0)) ' files in this directory \r '])
    fprintf('check the simulation options that are specified the options tab of the input spreadsheet \r')
    differentnumberoffiles = 1;
end

L = length(info0);
for i = 1:L
    s0 = info0(i).bytes;
    n0 = info0(i).name;
    for j = 1:length(info1)
        k = strcmp(info1(j).name, n0);
        if k, break, end
    end
    if k
        s1 = info1(j).bytes;
        if ~(s0==s1)
            fprintf(['\n Warning: the file size of ' info0(i).name ' is different from the verification output \r'])
            fprintf(['(' num2str(s1) ' instead of ' num2str(s0) ' bytes) \r'])
            differentsize = 1;
        end
        if (~strcmp(info0(i).name,'pars_and_input.csv') && ~strcmp(info0(i).name,'pars_and_input_short.csv'))
            D0 = dlmread([path0_ info0(i).name],',',2,0);
            D1 = dlmread([path1_ info1(j).name],',',2,0);
        elseif strcmp(info0(i).name,'pars_and_input_short.csv')
            continue
        else
            D0 = dlmread([path0_ info0(i).name],',',1,0);
            D1 = dlmread([path1_ info1(j).name],',',1,0);
        end
        if size(D0) == size(D1)
            if (sum(sum(D0-D1).^2))>1E-9
                fprintf(['\nWarning: data in the output file ' info0(i).name ' are different from the verification output \r '])
                h0 = textread([path0_ info0(i).name],'%s','bufsize', 1E9); %#ok<DTXTRD>
                spn = ceil(sqrt(size(D0,2)));
                figure(i)
                if spn>7
                    nr = min(size(D1, 1), size(D0, 1));
                    for z = 1:nr
                        plot(D0(z,:)','k'), hold on, plot(D1(z,:)','r')
                    end
                    title(info0(i).name, 'interpreter', 'none')
                else
                    h0 = strsplit(h0{1}, ',');
                    for m = 1:size(D0,2)
                        subplot(spn,spn,m)
                        plot(D0(:,m),'k'), hold on, plot(D1(:,m),'r')
                        title([info0(i).name h0(m)], 'interpreter', 'none')
                    end
                end
                differentcontent = 1;
            end
           % differentcontent = 1;
        end
        %end
    else
        fprintf(['\nWarning: the file ' info0(i).name ' was not found in the output\r'])
    end
end

if differentsize
    fprintf('\nWarning The size of some of the output files is different from the verification data \r')
    if differentcontent
        fprintf('Check if the startdate and enddate in the spreadsheet\r')
        fprintf('and the verification data in  are specified in "Dataset_Dir" in the Filenames tab of the input data spreadsheet \r')
    else
        fprintf('but the numerical values are the same. Possible cause: Different Matlab version\r')
    end
end
if ~(differentsize || differentcontent || differentnumberoffiles)
    fprintf('The output is the same as in the verification data set \r')
end
return