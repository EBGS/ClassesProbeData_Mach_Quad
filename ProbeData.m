classdef ProbeData < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here


    properties(SetAccess = protected, GetAccess = public)
        ShotNumber;
        ShotState;
        ProbeNumber;
        ElectronicBlockNumber;
        IndexPosition;      % Если используется много зондов одновремено, 
        % то следует каждому зонду приписать каждому зонду текущий индекс,
        % для извлечения о нем информацию из базы данных по положению Z и H
        
        ProbePositionAxisZ;
        ProbePositionAxisH;

        ProbeState;
        IsSignalsCorrected = 0;
        IsSignalsCalibrated = 0;
        PlasmaGunState;
        BeamState;

        BeamDelay;
        LimiterVoltage;
        ReciverVoltage;
        PlasmaGunVoltage;

        BeamDuration = 4;
        PlasmaGunDuration = 2.6;
        PlasmaGunDurationH2 = 2.8;
    end


    methods (Access = public)
        
        function obj = ProbeData(ShotNumber,ProbeNumber,ElectronicBlockNumber,IndexPosition)
            obj.ShotNumber = ShotNumber;
            obj.ProbeNumber = ProbeNumber;
            obj.IndexPosition = IndexPosition;
            obj.ElectronicBlockNumber = ElectronicBlockNumber;
        end

        function ExperimentalTable = GetExperimentalParameters(obj)
            cd TableExperimentData\
            load TableExperimentData.mat T
            cd ../


           ObjectNumber = max(size(obj));
           for i = 1:ObjectNumber
            index = (T.ShotNumber == obj(i).ShotNumber);
            IsShot = sum(index);

            if IsShot==1
                ShotNumber_(i) = obj(i).ShotNumber;
                BeamState_(i) = obj(i).BeamState;
                BeamDelay_(i) = obj(i).BeamDelay;
                LimiterVoltage_(i) = obj(i).LimiterVoltage;
                ReciverVoltage_(i) = obj(i).ReciverVoltage;
                PlasmaGunVoltage_(i) = obj(i).PlasmaGunVoltage;

                ProbePositionAxisZ_(i) = obj(i).ProbePositionAxisZ;
                ProbePositionAxisH_(i) = obj(i).ProbePositionAxisH;

                BeamDuration_(i) = obj(i).BeamDuration;
                PlasmaGunDuration_(i) = obj(i).PlasmaGunDuration;
                PlasmaGunDurationH2_(i) = obj(i).PlasmaGunDurationH2;
            else
                ExperimentalTable = 'No data entered in the table';
            end
           end

           ShotNumber_ = ShotNumber_';
           ProbePositionAxisZ_ = ProbePositionAxisZ_';
           ProbePositionAxisH_ = ProbePositionAxisH_';
           BeamState_ = BeamState_';
           BeamDelay_ = BeamDelay_';
           BeamDuration_ = BeamDuration_';
           LimiterVoltage_ = LimiterVoltage_';
           ReciverVoltage_ = ReciverVoltage_';
           PlasmaGunVoltage_ = PlasmaGunVoltage_';
           PlasmaGunDuration_ = PlasmaGunDuration_';
           PlasmaGunDurationH2_ = PlasmaGunDurationH2_';


           ExperimentalTable = table(ShotNumber_,ProbePositionAxisZ_,ProbePositionAxisH_,...
               BeamState_,BeamDelay_,BeamDuration_,LimiterVoltage_,ReciverVoltage_, ...
               PlasmaGunVoltage_,PlasmaGunDuration_,PlasmaGunDurationH2_);

        end

    end


    methods(Abstract = true)

        CorrectSignals(obj);
        DrawPlot(YLim1,YLim2);

    end


    methods (Access = protected)

        function obj = AddExperimentalParameters(obj)

            cd TableExperimentData\
            load TableExperimentData.mat T
            cd ../

            index = (T.ShotNumber == obj.ShotNumber);
            IsShot = sum(index);

            
            if IsShot==1
                T1 = T(index,1:end);

                obj.BeamState = T1.BeamState;
                obj.BeamDelay = T1.LimiterVoltage;
                obj.LimiterVoltage = T1.LimiterVoltage;
                obj.ReciverVoltage = T1.ReciverVoltage;
                obj.PlasmaGunVoltage = T1.PlasmaGunVoltage;

                if obj.IndexPosition == 1
                    obj.ProbePositionAxisZ = T1.ProbePositionAxisZ1;
                    obj.ProbePositionAxisH = T1.ProbePositionAxisH1;
                elseif obj.IndexPosition == 2
                    obj.ProbePositionAxisZ = T1.ProbePositionAxisZ2;
                    obj.ProbePositionAxisH = T1.ProbePositionAxisH2;
                end
            end
        end

    end


end