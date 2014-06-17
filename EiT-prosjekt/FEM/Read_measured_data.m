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
%  Loading data from validation script (THICK):   |
%-------------------------------------------------|
Data = load('../Data/validation_thick.mat');
forces = [0 200 400 600 800 1000];
disps1 = zeros(6,1);
for i=1:(length(disps1)-1)
    disps1(i+1) = 10^9*max(abs(Data.topdisp(i,:)));
end
clear Data;

%-------------------------------------------------|
%       Loading data from second validation:      |
%-------------------------------------------------|
Data = load('../Data/validation2.mat');
disps2 = zeros(6,1);
for i=1:(length(disps2)-1)
    disps2(i+1) = 10^9*max(abs(Data.topdisp(i,:)));
end
clear Data;



%-------------------------------------------------|
%       Getting out comparison plot:              |
%-------------------------------------------------|
figure
plot(disps1,forces,'ro-')
hold on;
plot(disps2,forces,'go-')
plot(beadA(:,2),beadA(:,3),'b')
xlabel('Displacement [nm]')
ylabel('Force [\muN]')
legend('Simulated (200nm shell)','Simulated (100nm shell)','Measured')
grid on;
axis([0 500 0 1000]);
