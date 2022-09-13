%"Создаем объекты"
clc;
clear all;
close all;
load mach_data.mat

%% "Определение распределений - Отношение токов - Ratio" 
close all;
t1 = 1.2; t2 = 2.5; RatioMin = 0; RatioMax = 50;

[mu,sigma,h] = GetCurrentRatio(mach_NB_6403_6413,[t1 t2],[RatioMin RatioMax]);
distr_ratio(1) = ProbeDataDistribution(h,mu,sigma,[6403 6413],[1.2 2.5]);
[mu,sigma,h] = GetCurrentRatio(mach_NB_6417_6427,[t1 t2],[RatioMin RatioMax]);
distr_ratio(2) = ProbeDataDistribution(h,mu,sigma,[6417 6427],[1.2 2.5]);

t1 = 2.85; t2 = 3.0; RatioMin = 0; RatioMax = 3;

[mu,sigma,h] = GetCurrentRatio(mach_NB_6403_6413,[t1 t2],[RatioMin RatioMax]);
distr_ratio(3) = ProbeDataDistribution(h,mu,sigma,[6403 6413],[1.2 2.5]);
[mu,sigma,h] = GetCurrentRatio(mach_NB_6417_6427,[t1 t2],[RatioMin RatioMax]);
distr_ratio(4) = ProbeDataDistribution(h,mu,sigma,[6417 6427],[1.2 2.5]);

%% "Определение распределений - токов I1, I2"
close all;
t1 = 1.2; t2 = 2.5; CurrentMin = 0; CurrentMax = 1;
[mu1, sigma1, mu2, sigma2, h] = GetCurrent(mach_NB_6403_6413,[t1 t2],[CurrentMin CurrentMax]);
distr_current_I1(1) = ProbeDataDistribution(h, mu1, sigma1,[6403 6413],[1.2 2.5]);
distr_current_I2(1) = ProbeDataDistribution(h, mu2, sigma2, [6403 6413],[1.2 2.5]);

[mu1, sigma1, mu2, sigma2, h] = GetCurrent(mach_NB_6417_6427,[t1 t2],[CurrentMin CurrentMax]);
distr_current_I1(2) = ProbeDataDistribution(h, mu1, sigma1, [6417 6427],[1.2 2.5]);
distr_current_I2(2) = ProbeDataDistribution(h, mu2, sigma2, [6417 6427],[1.2 2.5]);



t1 = 2.85; t2 = 3.0; CurrentMin = 0; CurrentMax = 0.1;
[mu1, sigma1, mu2, sigma2, h] = GetCurrent(mach_NB_6403_6413,[t1 t2],[CurrentMin CurrentMax]);
distr_current_I1(3) = ProbeDataDistribution(h, mu1, sigma1, [6403 6413],[1.2 2.5]);
distr_current_I2(3) = ProbeDataDistribution(h, mu2, sigma2, [6403 6413],[1.2 2.5]);

[mu1, sigma1, mu2, sigma2, h] = GetCurrent(mach_NB_6417_6427,[t1 t2],[CurrentMin CurrentMax]);
distr_current_I1(4) = ProbeDataDistribution(h, mu1, sigma1,[6417 6427],[1.2 2.5]);
distr_current_I2(4) = ProbeDataDistribution(h, mu2, sigma2,[6417 6427],[1.2 2.5]);







%% "Сохранение"
clearvars -except distr_current_I1 distr_current_I2 distr_ratio
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

load mach_data.mat
save('mach_data','distr_current_I1', 'distr_current_I2', 'distr_ratio','-append') 










