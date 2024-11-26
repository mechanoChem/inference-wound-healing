function lossGroup33rollingwin5F200000 = importLossFile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  LOSSGROUP33ROLLINGWIN5F200000 = IMPORTFILE(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the numeric
%  data.
%
%  LOSSGROUP33ROLLINGWIN5F200000 = IMPORTFILE(FILE, DATALINES) reads
%  data for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  lossGroup33rollingwin5F200000 = importfile("C:\Users\pkinn\Dropbox (University of Michigan)\vsiglobus\ht_migration\results\VSI_gamma_matrix\group2_A04_to_A06\Physics_Based_Time_Independent\loss_Group_3_3_rolling_win5_F200000.dat", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Dec-2022 09:46:47

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = [" ", ","];

% Specify column names and types
opts.VariableNames = ["e02", "e02_1", "e02_2", "e02_3", "e02_4", "e02_5", "e02_6", "e02_7"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
lossGroup33rollingwin5F200000 = readtable(filename, opts);

%% Convert to output type
lossGroup33rollingwin5F200000 = table2array(lossGroup33rollingwin5F200000);
end