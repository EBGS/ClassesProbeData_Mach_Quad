% Создаем таблицу данных - с предарительными данными
% Заполняем таблицу данными
% Очищаем идентичные строки
% Очищаем строки с одинаковыми выстрелами
% Извлекаем строку из таблицы с заданным нами номером выстрела

%% Генерируем таблицу данных по эксперименту
clear T;
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);

ShotNumber = [6618;6619];
BeamState = [1;1];
BeamDelay = [0;0];
LimiterVoltage = [150;150];
ReciverVoltage = [-200;-200];
PlasmaGunVoltage = [2.7;2.7];

ProbePositionAxisZ1 = [86;86];
ProbePositionAxisH1 = [95;95];
ProbePositionAxisZ2 = [0;0];
ProbePositionAxisH2 = [0;0];

T = table(ShotNumber,BeamState,BeamDelay,...
    LimiterVoltage,ReciverVoltage,PlasmaGunVoltage,...
    ProbePositionAxisZ1, ProbePositionAxisH1,ProbePositionAxisZ2,ProbePositionAxisH2);

cd TableExperimentData\
clearvars -except T
save TableExperimentData.mat T
clear T
cd ../
%% Добавляем в таблицу данные по эксперименту 
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);
cd TableExperimentData\
load TableExperimentData.mat T
cd ../

row = max(size(T.ShotNumber));

N = 6;
ShotNumber = [6307 6308 6312 6309 6310 6311]';
BeamState = ones(N,1)*1;
BeamDelay = ones(N,1)*0;
LimiterVoltage = ones(N,1)*(-150);
ReciverVoltage = ones(N,1)*(-200);
PlasmaGunVoltage = ones(N,1)*2.7;

ProbePositionAxisZ1 = ones(N,1)*(86);
ProbePositionAxisH1 = [0 0 0 0 0 0]';
ProbePositionAxisZ2 = ones(N,1)*0;
ProbePositionAxisH2 = ones(N,1)*0;

dT = table(ShotNumber,BeamState, BeamDelay,...
    LimiterVoltage,ReciverVoltage,PlasmaGunVoltage,...
    ProbePositionAxisZ1, ProbePositionAxisH1, ProbePositionAxisZ2,ProbePositionAxisH2);

T = [T;dT];

cd TableExperimentData\
clearvars -except T
save TableExperimentData.mat T
clear T
cd ../
%% Удаляем одинаковые копии
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);
cd TableExperimentData\
load TableExperimentData.mat T
cd ../

T = unique(T);

cd TableExperimentData\
clearvars -except T
save TableExperimentData.mat T
clear T
cd ../
%% Удаляем строки с одинаковыми номерами выстрелов
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);
cd TableExperimentData\
load TableExperimentData.mat T
cd ../

[newShotNumer,iT] = unique(T.ShotNumber);
newTable = T(iT,:);

cd TableExperimentData\
clearvars -except newTable
T = newTable;
save TableExperimentData.mat T
clear T newTable
cd ../
%% Сортируем таблицу
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);
cd TableExperimentData\
load TableExperimentData.mat T
cd ../

T = sortrows(T,1);

cd TableExperimentData\
clearvars -except T
save TableExperimentData.mat T
clear T
cd ../
%% Ищем нужный элемент
MainFolder = 'D:\Евгений\TheClassProbeData';
cd(MainFolder);
cd TableExperimentData\
load TableExperimentData.mat T
cd ../

ShotNumber = 6309;
index = (T.ShotNumber == ShotNumber);
IsShot = sum(index);

if IsShot==1
T1 = T(index,1:end)
T1.BeamState
T1.LimiterVoltage
end



%% Конец