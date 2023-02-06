function [Samples,Features,lgAUC,MutationMatrix] = parsexls2(file, sheet, col_lgAUC, col_Mut)
% col is the col-th data column containing lgAUC but 
% excluding the cell-line column 

%% read in cls file
[num,txt,raw]=xlsread(file, sheet);
%% Parse

% read name of cell-lines
Samples = txt(2:end,1);

% read log AUC values
lgAUC = num(:,col_lgAUC);

% read mutation matrix
MutationMatrix = num(:,col_Mut:end);

% read name of genetic features
Features = txt(1,(col_Mut+1):end);