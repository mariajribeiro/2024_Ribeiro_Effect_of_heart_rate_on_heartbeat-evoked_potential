x = importdata('check-Rpeak_time_sub-AB58_task-passive_run1.mat');
y = importdata('Tpeak_time_sub-AB58_task-passive_run1.mat');
fileID = fopen('check-Rpeak_Tpeak_sub-AB58_task-passive_run1.txt','wt');
fprintf(fileID,'%s %s\n','latency','type');
for i = 1:length(x)
fprintf(fileID,'%f %s\n', x(i),'Rpeak');
end
for i = 1:length(y)
fprintf(fileID,'%f %s\n', y(i),'Tpeak');
end
fclose(fileID);

%only R peak events:AB21,50,55
x = importdata('check-Rpeak_time_sub-AB55_task-passive_run1.mat')
fileID = fopen('check-Rpeak_sub-AB55_task-passive_run1.txt','wt');
fprintf(fileID,'%s %s\n','latency','type');
for i = 1:length(x)
fprintf(fileID,'%f %s\n', x(i),'Rpeak');
end
fclose(fileID);