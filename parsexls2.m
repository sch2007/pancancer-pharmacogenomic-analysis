function [Samples,Features,AUC,MutationMatrix] = parsexls2(file, sheet, col_AUC, col_Mut)
% col is the col-th data column containing AUC but 
% excluding the cell-line column 

%% read in cls file
[num,txt,raw]=xlsread(file, sheet);
%% Parse

% read name of cell-lines
Samples = txt(2:end,1);

% read log AUC values
AUC = num(:,col_AUC);

% read mutation matrix
MutationMatrix = num(:,col_Mut:end);

% read name of genetic features
Features = txt(1,(col_Mut+1):end);