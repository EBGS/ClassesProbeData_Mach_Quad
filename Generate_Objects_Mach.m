%% "Создаем файл для сохранения объектов mach_data"
clc;
clear all;
close all;
MainFolder = 'D:\Евгений\TheClassProbeData';
save('mach_data','MainFolder') % "Только при первой инициализации сохранения"
%%
clc;
clear all;
close all;
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

% "NB6403-6407,6409-6413"
ShotNumber = [6403 6404 6405 6406 6407 6409 6410 6411 6412 6413];
num = max(size(ShotNumber));
for i=1:num
mach_NB_6403_6413(i) = ProbeDataMach(ShotNumber(i),1,1,1);
end
mach_NB_6403_6413.SetGoal('Профиль зонда Маха без нейтрального пучка z = -136 см')



% "NB6417-6427"
ShotNumber = [6403 6404 6405 6406 6407 6409 6410 6411 6412 6413];
num = max(size(ShotNumber));
for i=1:num
mach_NB_6417_6427(i) = ProbeDataMach(ShotNumber(i),1,1,1);
end
mach_NB_6417_6427.SetGoal('Профиль зонда Маха с нейтрльным пучком z = -136 см')

%% "Чистим сигналы - убираем эффекты блока Электроники"
mach_NB_6403_6413.CorrectSignals();
mach_NB_6417_6427.CorrectSignals();
%% Калибруем сигналы
lambda = 0.7497; % sigma = 0.1376;
mach_NB_6403_6413.SetCalibration(lambda)
mach_NB_6417_6427.SetCalibration(lambda)







%% "Рисуем графики"
close all;
mach_NB_6403_6413.DrawPlot([-0.5 4],[0 0.75])
mach_NB_6417_6427.DrawPlot([-0.5 4],[0 0.75])
%% "Рисуем отношение токов"
close all;
mach_NB_6403_6413.DrawPlotRatio([-0.5 4],[0 25])
mach_NB_6417_6427.DrawPlotRatio([-0.5 4],[0 25])
%% "Рисуем отношение токов - во временном интервале"
close all;
t1 = 2.85; t2 = 3.0; RatioMin = 0; RatioMax = 3;
DrawPlotRatioInterval(mach_NB_6403_6413,[t1 t2],[RatioMin RatioMax],"pdf");
DrawPlotRatioInterval(mach_NB_6417_6427,[t1 t2],[RatioMin RatioMax],"pdf");






%% "Сохранение"
clearvars -except mach_NB_6403_6413 mach_NB_6417_6427
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

load mach_data.mat
save('mach_data','mach_NB_6403_6413','mach_NB_6417_6427','-append') 
















