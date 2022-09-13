%"Создаем объекты"
% clc;
% clear all;
% close all;
% load quad_data.mat

%% "Определение распределений - Электрическое поле" 
close all;
VoltageShift = -6;
[h,Er,dEr,t1,t2] = quad_NB_6097_6112.Get_VAC_Er(5,VoltageShift);
distr_Er(1) = ProbeDataDistribution(h, Er, dEr,[t1 t2],[6403 6413]);
[h,Er,dEr,t1,t2] = quad_NB_6097_6112.Get_VAC_Er(6,VoltageShift);
distr_Er(2) = ProbeDataDistribution(h, Er, dEr,[t1 t2],[6403 6413]);
[h,Er,dEr,t1,t2] = quad_NB_6097_6112.Get_VAC_Er(7,VoltageShift);
distr_Er(3) = ProbeDataDistribution(h, Er, dEr,[t1 t2],[6403 6413]);

%% "Отзеркаливание профилей"
GetAntiSymmetricDistribution(distr_Er,[155 180]);
distr_Er.DrawDistributionFit([-50 50],[-30 30],0.75,'Er','h, мм','Er, В',1)

%%
% %% "Сохранение"
% clearvars -except distr_current_I1 distr_current_I2 distr_ratio
% MainFolder = 'D:\Евгений\The Class ProbeData';
% cd(MainFolder);
% 
% load mach_data.mat
% save('mach_data','distr_current_I1', 'distr_current_I2', 'distr_ratio','-append') 
% save('mach_data','distr_current_I1', 'distr_current_I2', 'distr_ratio') 
% 
% 







