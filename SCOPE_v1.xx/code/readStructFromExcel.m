function data = readStructFromExcel(filename, sheetName, headerIdx, dataIdx, data_is_char, data_in_rows)
% Read data into a struct with names matching those found in the first column/row
%  default is for data to be in columns (A and B); if data_in_rows = true, data are in rows 1 & 2
%  example:
%   readStructFromExcel('../input_data.xlsx', 'options', 3, 1)
%   readStructFromExcel('../input_data.xlsx', 'filenames', 1, 2, true)

% default values:
if nargin < 3
    headerIdx = 1;
end
if nargin < 4
    dataIdx = 2;
end
if nargin < 5
    data_is_char = false;
end
if nargin < 6
    data_in_rows = false;
end

% read in the spreadsheet
% general note: work with all_as_cell to keep rows & columns in sync; MATLAB does NOT keep data and texts aligned
%   (readtable only slightly better - it works but will convert mixed columns to string)
if data_in_rows
    %NOTE: 'basic' is compatible with Mac but MATLAB (2013b) complains it can't read Unicode '.xls' in basic mode
    %  Solution: save as xlsx ?!
    [~, ~, allCells] = xlsread(filename, sheetName, '');  % = [data, texts, allCells]
else
    [~, ~, allCells] = xlsread(filename, sheetName, '');  % = [data, texts, allCells]
    allCells = allCells'; % transpose so data are in the same column as headers
end
%  data are now in columns

% delete empty columns
 % define two "helper functions" for eliminating null entries
isNotNumeric = @(x) cellfun(@(y) ischar(y) | any(isnan(y)), x);  % any is needed because matlab treats string as char array
isCharCell = @(x) cellfun(@(y) ischar(y), x);  

validHeaders = arrayfun(isCharCell, allCells(headerIdx, :));
if data_is_char
    validData = arrayfun(isCharCell, allCells(dataIdx, :));  % , 'UniformOutput', true
else % numeric data:
    validData = ~arrayfun(isNotNumeric, allCells(dataIdx, :));  % , 'UniformOutput', true
end

dataCells = allCells([headerIdx, dataIdx], validHeaders & validData);

for idx = 1:size(dataCells,2)
    varName = strrep(dataCells{1, idx}, ' ', '');
    data.(varName) = dataCells{2, idx};
end