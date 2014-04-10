%--------------------------------------------------|
%       Loading data from Excel spreadsheet:       |
%--------------------------------------------------|

Data = xlsread('../Data/De 4.1.3 cyclic/Plots.xlsx');
%Column 1: Time (s)
%Column 2: Displacement (nm),
%Column 3: Force (\mu N)
% Cols = [1 5 9 14 19 24];
% Beads = ['A' 'B' 'C' 'D' 'E' 'F'];
beadA = getBeadData(Data,1);
beadB = getBeadData(Data,5);
beadC = getBeadData(Data,9);
beadD = getBeadData(Data,14);
beadE = getBeadData(Data,19);
beadF = getBeadData(Data,24);
clear Data;

%-------------------------------------------------|
%       Loading data from validation script:      |
%-------------------------------------------------|
Data = load('../Data/validation.mat');
forces = [0 200 400 600 800 1000];
disps = zeros(6,1);
for i=1:(length(disps)-1)
    disps(i+1) = 10^9*max(abs(Data.topdisp(i,:)));
end
clear Data;

%-------------------------------------------------|
%       Getting out comparison plot:              |
%-------------------------------------------------|
figure
plot(disps,forces,'ro-')
hold on;
plot(beadA(:,2),beadA(:,3),'b')
xlabel('Displacement [nm]')
ylabel('Force [\muN]')
legend('Simulated','Measured')
grid on;
axis([disps(1) max(beadA(:,2)) forces(1) forces(end)]);
