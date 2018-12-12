path = "D:\projects Python\SCOPE\SCOPE_v1.70\docs\get_input_structs\input_data - Copy.xlsx";
tab = readtable(path, 'Sheet', 'inputdata', 'Range', '6:90');
% tab = removevars(tab,{'Unit','Description'});
% delete text field (more generaly might be a check on cell arrays of fields
tab.Unit = [];
tab.Description = [];
% strcmp(tab.Input_for_SCOPE, 'tts')  % this works
% tab(strcmp(tab.Input_for_SCOPE, 'tts'))  % this does not work
tab(strcmp(tab.Input_for_SCOPE, 'tts'), 1)
tab(strcmp(tab.Input_for_SCOPE, 'tts'), :)  % complete row
val = tab(strcmp(tab.Input_for_SCOPE, 'tts'), 2);
% tab(strcmp(tab.Input_for_SCOPE, 'tts'), 2) + 10  % this is a table => not possible 
val_cell = tab{strcmp(tab.Input_for_SCOPE, 'tts'), 2};
val_char = tab{strcmp(tab.Input_for_SCOPE, 'tts'), 2}{:};
vals = tab{strcmp(tab.Input_for_SCOPE, 'tts'), 3:7};  % Var2, Var7 cells of chars
val = str2num(val_char);

%% simple
% read skipping rows to have headers
tab = readtable(path, 'Sheet', 'inputdata', 'Range', '6:90');

% delete text field (more generaly might be a check on cell arrays of fields
tab.Unit = [];
tab.Description = [];
tab.Var10 = [];

% get values
values = tab{strcmp(tab.Variable, 'tts'), 2:end};
