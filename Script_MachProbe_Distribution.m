
%% "Создаем объекты"
clc;
clear all;
close all;
load mach_data.mat

%% "Отзеркаливание профилей"
GetSymmetricDistribution(distr_ratio,[55 80]);
%%
close all; figNumber1 = 1; figNumber2 = 3; 
distr_ratio(logical([1 1 0 0])).DrawDistributionFit([-50 50],[0 25],0.6,'Отношение токов 1 и 2','h, мм', 'ratio',figNumber1)
distr_ratio(logical([0 0 1 1])).DrawDistributionFit([-50 50],[0 2],0.6,'Отношение токов 1 и 2','h, мм', 'ratio',figNumber2)
%%
close all;
distr_ratio(logical([1 1 0 0])).DrawDistribution([45 105],[0 50],1,'Отношение токов 1 и 2','h, мм', 'ratio',1)
distr_ratio(logical([0 0 1 1])).DrawDistribution([45 105],[0 5],2,'Отношение токов 1 и 2','h, мм', 'ratio',1)
distr_ratio(logical([1 1 1 1])).DrawDistribution([45 105],[0 50],3,'Отношение токов 1 и 2','h, мм', 'ratio',1)
%%
close all;
distr_ratio(logical([1 1 0 0])).DrawErrorBarDistribution([45 105],[0 40],1,'Отношение токов 1 и 2','h, мм', 'ratio')
distr_ratio(logical([0 0 1 1])).DrawErrorBarDistribution([45 105],[0 3],2,'Отношение токов 1 и 2','h, мм', 'ratio')
distr_ratio(logical([1 1 1 1])).DrawErrorBarDistribution([45 105],[0 50],3,'Отношение токов 1 и 2','h, мм', 'ratio')

% %%
% obj.SetGoal('Профиль зонда Маха без пучка')
% %%
% goal1 = GetGoal(obj)
%%
[mu1, sigma1, mu2, sigma2, h] = GetCurrentDistribution(obj,[1.2 2.5],[-1 1]);
close all;
distr(1) = ProbeDataDistribution(h,mu1,sigma1,[6403 6413],[1.2 2.5]);
distr(2) = ProbeDataDistribution(h,mu2,sigma2,[6417 6427],[1.2 2.5]);
%%
distr(logical([1 1])).DrawDistribution([45 105],[0 5],1,'Токи 1 и 2','h, мм', 'I, А',1)








