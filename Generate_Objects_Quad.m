% "Создаем файл для сохранения объектов mach_data"
clc;
clear all;
close all;
MainFolder = 'D:\Евгений\TheClassProbeData';
save('quad_data','MainFolder') % "Только при первой инициализации сохранения"
%% Создаем объекты
clc;
clear all;
close all;
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

% "NB6118-6127,6133-6137"
ShotNumber = [6118 6119 6120 6121 6122 6123 6124 6125 6126 6127 6133 6134 6135 6136 6137];
num = max(size(ShotNumber));
for i=1:num
quad_NB_6117_6137(i) = ProbeDataQuad(ShotNumber(i),1,1,1);
end
quad_NB_6117_6137.SetGoal('Профиль Четырех-электродного зонда без нейтрального пучка z = 86 см')

% "NB6097-6098,6100-6112"
ShotNumber = [6097 6098 6100 6101 6102 6103 6104 6105 6106 6107 6108 6109 6110 6111 6112];
num = max(size(ShotNumber));
for i=1:num
quad_NB_6097_6112(i) = ProbeDataQuad(ShotNumber(i),1,1,1);
end
quad_NB_6097_6112.SetGoal('Профиль Четырех-электродного зонда с пучком задержка tau = 0 мс z = 86 см')


% "NB6139-6140,6142-6147,6152,6154-6158,6160"
ShotNumber = [6139 6140 6142 6143 6144 6145 6146 6147 6152 6154 6155 6156 6157 6158 6160];
num = max(size(ShotNumber));
for i=1:num
quad_NB_6139_6160(i) = ProbeDataQuad(ShotNumber(i),1,1,1);
end
quad_NB_6139_6160.SetGoal('Профиль Четырех-электродного зонда с пучком задержка tau = 0.8 мс z = 86 см')


%% "Чистим сигналы - убираем эффекты блока Электроники"
quad_NB_6117_6137.CorrectSignals();
quad_NB_6097_6112.CorrectSignals();
quad_NB_6139_6160.CorrectSignals();

%% Делаем автокалибровку сигналов
quad_NB_6117_6137.SetAutoCalibration();
quad_NB_6097_6112.SetAutoCalibration();
quad_NB_6139_6160.SetAutoCalibration();


% % % % % %% Рисуем Нормированный ионный ток (без калибровки)
% % % % % close all;
% % % % % quad_NB_6117_6137(2).DrawCurrentNormilized([0.0 2])
% % % % % %% Величина Нормированного ионного тока (только для одиночных объектов)
% % % % % close all;
% % % % % [mu, sigma,time] = quad_NB_6097_6112(2).GetCurrentNormilizedIonParthTimePoints();
% % % % % %% Величина Нормированного ионного тока (средняя в выстреле в разных выстрелах)
% % % % % close all;
% % % % % [mu, sigma] = quad_NB_6097_6112.GetCurrentNormilizedIonParth();


%% Получаем все ВАХ во всех выстрелах Uexp{s,j} s - номер выстела, j - номер временного интервала
[Uexp,Jexp,t1,t2,tm,Er,Isat] = quad_NB_6097_6112.GetVac();
plot(Uexp{1,12},Jexp{1,12},'.')
ylim([-2 2])

%% Получаем параметры плазмы Te, n, U3, Er
[Uexp_reduce,Jexp_reduce,X_Fit,Y_Fit,n,Te,U3,dn,dTe,dU3] = GetPlasmaParameters(quad_NB_6097_6112,[-28 7],[-5 1.5]);
%% Рисуем ВАХ
close all
DrawVAC(quad_NB_6097_6112(5),[-30 10],[-3 2],[1 14])
%%
plot(Uexp_reduce{5,3},Jexp_reduce{5,3},'.',X_Fit{5,3},Y_Fit{5,3},'-')
%% "Рисуем графики"
% close all;
% %quad_NB_6117_6137.DrawPlot([-0.05 0.15],[-50 50])
% %quad_NB_6097_6112.DrawPlot([-0.05 0.15],[-50 50])
% quad_NB_6139_6160.DrawPlot([-0.03 0.15],[-50 50])
% %% Таблица Экспериментальных данных
table = quad_NB_6117_6137(7).GetExperimentalParameters()
%% "Сохранение"
clearvars -except quad_NB_6117_6137 quad_NB_6097_6112 quad_NB_6139_6160
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

% load quad_data.mat
% save('quad_data','quad_NB_6117_6137','quad_NB_6097_6112','quad_NB_6139_6160','-append') 
save('quad_data','quad_NB_6117_6137','quad_NB_6097_6112','quad_NB_6139_6160') 


