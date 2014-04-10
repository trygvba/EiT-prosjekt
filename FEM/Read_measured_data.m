Data = xlsread('../Measurements_results/Plots.xlsx');

%Column 1: Time (s)
%Column 2: Displacement (nm),
%Column 3: Force (\mu N)
beadA = Data(find(~isnan(Data(:,1))),1:3);
beadB = Data(find(~isnan(Data(:,5))),5:7);