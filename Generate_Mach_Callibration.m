%% "Отбраковываем выстрелы"
% исходные выстрелы [6307 6308 6312]
% исходные выстрелы [6309 6310 6311]
%% "Отбраковываем выстрелы"
% исходные выстрелы [6307 6308 6312]
% исходные выстрелы [- 6310 -]
%% "Создаем объекты для Калибровки"
clc;
clear all;
close all;
cal1(1) = ProbeDataMach(6307,1,1,1);
cal1(2) = ProbeDataMach(6308,1,1,1);
cal1(3) = ProbeDataMach(6312,1,1,1);
cal2(1) = ProbeDataMach(6310,1,1,1);
cal2(2) = ProbeDataMach(6311,1,1,1);
cal2(3) = ProbeDataMach(6309,1,1,1);
%% "Корректируем сигналы"
cal1.CorrectSignals();
cal2.CorrectSignals();
%% "Рисуем графики"
close all;
cal1.DrawPlot([0 4.2],[0 0.2])
cal2.DrawPlot([0 4.2],[0 0.2])
%% "Получаем величины отношений токов - ratio"
close all;
[mu1 , sigma1] = GetCurrentRatio(cal1,[1.2 2.5],[0 2.5]);
[mu2 , sigma2] = GetCurrentRatio(cal2,[1.2 2.5],[0 2.5]);
%% "Получаем четверное отношений токов - lambda"
close all;
[lambda, sigma] = GetCallibrationLambdaValue(cal1,cal2,[1.2 2.5],[0 2]);
%% "Установить калибровку"
cal1.SetCalibration(lambda)
cal2.SetCalibration(lambda)
%% "Рисуем графики"
close all;
cal1.DrawPlot([0 4.2],[0 0.2])
cal2.DrawPlot([0 4.2],[0 0.2])
%% "Получаем четверное отношений токов - lambda  (повторно) "
close all;
[lambda, sigma] = GetCallibrationLambdaValue(cal1,cal2,[1.2 2.5],[0 2]);








