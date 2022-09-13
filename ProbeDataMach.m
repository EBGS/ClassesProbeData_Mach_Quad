classdef ProbeDataMach < ProbeData & handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = protected, GetAccess = public)
        t;
        Im1;
        Im2;

        Lambda;

        tgun;
        Igun;

        Ratio;

        Goal;

    end


    methods (Access = public)

        function obj = ProbeDataMach(ShotNumber,ProbeNumber,ElectronicBlockNumber,IndexPosition)
            
            obj = obj@ProbeData(ShotNumber,ProbeNumber,ElectronicBlockNumber,IndexPosition);
            
            if (nargin>0)   
                % Определяем состояния эксеримента:
                obj = WhatShotState(obj);
                % ProbeState = 0 or 1
                cd Shots\
                if obj.ProbeState > 0
                    Im1_=importdata(['Im2.',num2str(ShotNumber)]);
                    Im2_=importdata(['Im1.',num2str(ShotNumber)]);
                    obj.Im1 = -(Im1_(1:end,2));
                    obj.Im2 = -(Im2_(1:end,2));
                    obj.t = (Im1_(1:end,1)/1000);
                    obj.Ratio = obj.Im2./obj.Im1;
                end

                if obj.PlasmaGunState > 0
                    Igun_=importdata(['Gun_Cath_Curr.',num2str(ShotNumber)]);
                    obj.Igun = (Igun_(1:end,2));
                    obj.tgun = (Igun_(1:end,1)/1000);
                end
            end
            cd ../

            if (nargin>0)
            obj = AddExperimentalParameters(obj);   
            end
          
        end

        function CorrectSignals(obj)
            ObjectNumer = max(size(obj));
            for s = 1:ObjectNumer
            if (obj(s).ProbeState == 1)&&(obj(s).IsSignalsCorrected == 0)
                % Коррекция временных сигналов
                obj(s).t = obj(s).t - 0.1;
                % Вычитание фона
                time = obj(s).t;
                Jm1 = obj(s).Im1;
                Jm2 = obj(s).Im2;
                Im1_ground=mean(Jm1(time>-0.2&time<-0.1));
                Im2_ground=mean(Jm2(time>-0.2&time<-0.1));
                obj(s).Im1=(obj(s).Im1-Im1_ground);
                obj(s).Im2=(obj(s).Im2-Im2_ground);
                obj(s).Ratio = obj(s).Im2./obj(s).Im1;
            end
            if obj(s).PlasmaGunState == 1
                % Коррекция временных сигналов
                obj(s).tgun = obj(s).tgun - 0.035;
                % Вычитание фона
                timeGun = obj(s).tgun;
                Jgun = obj(s).Igun;
                Igun_ground=mean(Jgun(timeGun>-0.2&timeGun<-0.1));
                obj(s).Igun=(obj(s).Igun-Igun_ground);              
            end
            obj(s).IsSignalsCorrected = 1;
            end
        end

        function SetCalibration(obj,Lambda)
            % задание коэффициентов коррекции Амплитуд
            K1=(Lambda).^(1/4);
            K2=(Lambda).^(-1/4);

            ObjectNumer = max(size(obj));
            for s = 1:ObjectNumer
                 if (obj(s).ProbeState == 1)&&(obj(s).IsSignalsCalibrated == 0)
                    % Коррекция Амплитуд сигналов
                    obj(s).Lambda = Lambda;
                    obj(s).Im1 = K1*obj(s).Im1;
                    obj(s).Im2 = K2*obj(s).Im2;
                    obj(s).Ratio = obj(s).Im2./obj(s).Im1;
                    obj(s).IsSignalsCalibrated = 1;
                 end
            end

        end

        function SetGoal(obj,Goal)
            ObjectNumer = max(size(obj));
            for s = 1:ObjectNumer
                obj(s).Goal = Goal;
            end
        end

        function goal = GetGoal(obj)
            goal = obj.Goal;
        end

        function DrawPlot(obj,limit_X,limit_Y)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
                if obj(s).ProbeState == 1
                    f = figure(obj(s).ShotNumber);
                    f.Color = [1 1 1];
                    f.Position = [170 50 1200 730];
                    p1 = plot(obj(s).t,obj(s).Im1,obj(s).t,obj(s).Im2);
                    p1(1).Color  = [1 0 0];
                    p1(2).Color  = [0 0 0];
                    p1(1).LineWidth = 1;
                    p1(2).LineWidth = 1;
                    s1  = gca;
                    s1.XLim = limit_X;
                    s1.YLim = limit_Y;
                    % s1.XGrid = 'on';
                    % s1.YGrid = 'on';
                    s1.Title.String = 'Токи на зондовые проволоки 1 и 2';
                    s1.XLabel.String = 't, мс';
                    s1.YLabel.String = 'I, А';
                    s1.XLabel.FontWeight = 'bold';
                    s1.YLabel.FontWeight = 'bold';
                    s1.FontSize = 16;
                    s1.FontName = 'Times New Roman';
                    s1.FontAngle = 'italic';
                    % s1.Color = [1 1 0.5];
                    legend('I_{1}','I_{2}')
                else
                    disp('There is no Shot in folder: D:\Евгений\The Class ProbeData\Shots');
                end
            end
        end

        function DrawPlotRatio(obj,limit_X,limit_Y)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
                if obj(s).ProbeState == 1

                    f = figure(obj(s).ShotNumber);
                    f.Color = [1 1 1];
                    f.Position = [170 50 1200 730];
                    p2 = plot(obj(s).t,obj(s).Ratio);
                    p2.Color  = [0 0 0];
                    p2.LineWidth = 1;
                    s2  = gca;
                    s2.XLim = limit_X;
                    s2.YLim = limit_Y;
                    % s2.XGrid = 'on';
                    % s2.YGrid = 'on';
                    s2.Title.String = 'Отношение токов на проволочки 1 и 2 в зонде Маха';
                    s2.XLabel.String = 't, мс';
                    s2.YLabel.String = 'I_{1} / I_{2}';
                    s2.XLabel.FontWeight = 'bold';
                    s2.YLabel.FontWeight = 'bold';
                    s2.FontSize = 16;
                    s2.FontName = 'Times New Roman';
                    s2.FontAngle = 'italic';
                    % s2.Color = [1 1 0.5];
                else
                    disp('There is no Shot in folder: D:\Евгений\The Class ProbeData\Shots');
                end
            end
        end

        function [lambda, sigma] = GetCallibrationLambdaValue(obj1,obj2,TimeInterval,RatioInterval)

            [lambda, sigma, time, Ratio1, Ratio2, Lambda, DataTimeInterval,DataLambdaInterval, x_exper,y_exper_cdf,x,y_cdf,y_pdf] = CalculateLambda(obj1,obj2,TimeInterval,RatioInterval);

            f1 = figure(1);
            f1.Color = [1 1 1];
            plot(x_exper,y_exper_cdf,'.')
            hold on;
            plot(x,y_cdf,'LineWidth',1.5)
            grid on
            title(['CDF of normal distribution with \mu = ',num2str(lambda),' and \sigma = ',num2str(sigma)])
            xlabel('\lambda')
            ylabel('CDF')
            ylim([0 1.05])

            f2 = figure(2);
            f2.Color = [1 1 1];
            nbins = 55;
            histogram(x_exper,nbins,'Normalization','pdf')
            hold on
            plot(x,y_pdf,'LineWidth',1.5)
            grid on
            title(['PDF of normal distribution with \mu = ',num2str(lambda),' and \sigma = ',num2str(sigma)])
            xlabel('\lambda')
            ylabel('PDF')
        end

        function DrawPlotRatioInterval(obj,TimeInterval,RatioInterval,DistributionType)
    
            ObjectNumber = max(size(obj));
            mu(ObjectNumber) = 0;
            sigma(ObjectNumber) = 0;
            for s = 1:ObjectNumber

                if (obj(s).ProbeState == 1)&&(DistributionType == 'pdf')
                    [mu(s), sigma(s),time,Ratio, DataTimeInterval, DataRatioInterval,x_exper,y_exper_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj(s),TimeInterval,RatioInterval);
                    f2 = figure(obj(s).ShotNumber);
                    f2.Color = [1 1 1];
                    nbins = 55;
                    histogram(x_exper,nbins,'Normalization','pdf')
                    hold on
                    plot(x,y_pdf,'LineWidth',1.5)
                    grid on
                    title(['PDF of normal distribution with \mu = ',num2str(mu(s)),' and \sigma = ',num2str(sigma(s))])
                    xlabel('\lambda')
                    ylabel('PDF')
                elseif (obj(s).ProbeState == 1)&&(DistributionType == 'cdf')
                    [mu(s), sigma(s),time,Ratio, DataTimeInterval, DataRatioInterval,x_exper,y_exper_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj(s),TimeInterval,RatioInterval);
                    f1 = figure(obj(s).ShotNumber);
                    f1.Color = [1 1 1];
                    plot(x_exper,y_exper_cdf,'.')
                    hold on;
                    plot(x,y_cdf,'LineWidth',1.5)
                    grid on
                    title(['CDF of normal distribution with \mu = ',num2str(mu(s)),' and \sigma = ',num2str(sigma(s))])
                    xlabel('\lambda')
                    ylabel('CDF')
                    ylim([0 1.05])
                else
                    disp('There is no Shot in folder: D:\Евгений\The Class ProbeData\Shots');
                end
            end

        end
   
        function [mu , sigma, h] = GetCurrentRatio(obj,TimeInterval,RatioInterval)
  
            ObjectNumber = max(size(obj));
            mu = zeros(1,ObjectNumber);
            sigma = zeros(1,ObjectNumber);
            h = zeros(1,ObjectNumber);
            
            for s = 1:ObjectNumber
                if obj(s).ProbeState == 1
                    [mu(s), sigma(s)] = CalculateCurrentRatio(obj(s),TimeInterval,RatioInterval);
                    h(s) = obj(s).ProbePositionAxisH;
                end
            end

        end

        function [mu1, sigma1, mu2, sigma2, h] = GetCurrent(obj,TimeInterval,RatioInterval)

            ObjectNumber = max(size(obj));
            mu1 = zeros(1,ObjectNumber);
            sigma1 = zeros(1,ObjectNumber);
            mu2 = zeros(1,ObjectNumber);
            sigma2 = zeros(1,ObjectNumber);
            h = zeros(1,ObjectNumber);

            for s = 1:ObjectNumber
                if obj(s).ProbeState == 1
                    [mu1(s), sigma1(s), mu2(s), sigma2(s)] = CalculateCurrent(obj(s),TimeInterval,RatioInterval);
                    h(s) = obj(s).ProbePositionAxisH;
                end
            end

        end

        
    end

    methods (Access = private)

        function obj = WhatShotState(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cd Shots\

            fileID_Im1 = fopen(['Im1.' num2str(obj.ShotNumber)]);
            fileID_Im2 = fopen(['Im2.' num2str(obj.ShotNumber)]);
            fileID_Igun = fopen(['Gun_Cath_Curr.',num2str(obj.ShotNumber)]);

            fileID_Probe = (fileID_Im1 >= 0)&(fileID_Im2 >= 0);
            fileID_Gun = (fileID_Igun >= 0);

            obj.ProbeState = fileID_Probe;
            obj.PlasmaGunState = fileID_Gun;
            obj.ShotState = (fileID_Gun&fileID_Probe);
            cd ../
        end

        function [mu1, sigma1, mu2, sigma2] = CalculateCurrent(obj,TimeInterval,RatioInterval)
            ObjectNumber = max(size(obj));
            time = obj(1).t;

            current1 = 0;
            current2 = 0;
            for i=1:ObjectNumber
                current1 = current1 + obj(i).Im1;
                current2 = current2 + obj(i).Im2;
            end
            Current1=current1/ObjectNumber;
            Current2=current2/ObjectNumber;

            t1 = TimeInterval(1);           t2 = TimeInterval(2);           % Интервал времени
            CurrentMin = RatioInterval(1);    CurrentMax = RatioInterval(2);    % Интервал отнношений
            DataTimeInterval = time(time>t1 & time<t2);
            DataCurrentInterval_1_ = Current1(time>t1 & time<t2);
            DataCurrentInterval_2_ = Current2(time>t1 & time<t2);


            DataCurrentInterval_1 = sort(DataCurrentInterval_1_);
            DataCurrentInterval_2 = sort(DataCurrentInterval_2_);
            DataCurrentInterval_1_Sorted = DataCurrentInterval_1(DataCurrentInterval_1>CurrentMin&DataCurrentInterval_1<CurrentMax);
            DataCurrentInterval_2_Sorted = DataCurrentInterval_2(DataCurrentInterval_2>CurrentMin&DataCurrentInterval_2<CurrentMax);

            N1 = max(size(DataCurrentInterval_1_Sorted));
            N2 = max(size(DataCurrentInterval_2_Sorted));


            mu1 = mean(DataCurrentInterval_1_Sorted);
            sigma1 = std(DataCurrentInterval_1_Sorted);

            mu2 = mean(DataCurrentInterval_2_Sorted);
            sigma2 = std(DataCurrentInterval_2_Sorted);

        end
        
        function [mu, sigma,time,Ratio, DataTimeInterval, DataRatioInterval,x_exper,y_exper_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj,TimeInterval,RatioInterval)
            ObjectNumber = max(size(obj));
            time = obj(1).t;

            ratio = 0;
            for i=1:ObjectNumber
                ratio = ratio + obj(i).Ratio;
            end
            Ratio=ratio/ObjectNumber;

            t1 = TimeInterval(1);           t2 = TimeInterval(2);           % Интервал времени
            RatioMin = RatioInterval(1);    RatioMax = RatioInterval(2);    % Интервал отнношений
            DataTimeInterval = time(time>t1 & time<t2);   DataRatioInterval = Ratio(time>t1 & time<t2);


            SortedRatio_ = sort(DataRatioInterval);
            SortedRatio = SortedRatio_(SortedRatio_>RatioMin&SortedRatio_<RatioMax);
            N = max(size(SortedRatio));

            x_exper = SortedRatio;
            y_exper_cdf = ((1:1:N)/N)';
            mu = mean(SortedRatio);
            sigma = std(SortedRatio);

            x = RatioMin:0.005:RatioMax;
            y_cdf = cdf('Normal',x,mu,sigma)';
            y_pdf = pdf('Normal',x,mu,sigma)';


        end

        function [lambda, sigma, time, Ratio1, Ratio2, Lambda, DataTimeInterval,DataLambdaInterval, x_exper,y_exper_cdf,x,y_cdf,y_pdf] = CalculateLambda(obj1,obj2,TimeInterval,RatioInterval)
            % Первый объект
            ObjectNumber1 = max(size(obj1));
            time = obj1(1).t;
            ratio1 = 0;
            for i=1:ObjectNumber1
                ratio1 = ratio1 + obj1(i).Ratio;
            end
            Ratio1 = ratio1/ObjectNumber1;
            % Второй объект
            ObjectNumber2 = max(size(obj2));
            ratio2 = 0;
            for i=1:ObjectNumber2
                ratio2 = ratio2 + obj2(i).Ratio;
            end
            Ratio2=ratio2/ObjectNumber2;
            % Lambda
            Lambda = Ratio1.*Ratio2;
            t1 = TimeInterval(1);           t2 = TimeInterval(2);           % Интервал времени
            RatioMin = RatioInterval(1);    RatioMax = RatioInterval(2);    % Интервал отнношений
            DataTimeInterval = time(time>t1 & time<t2);   DataLambdaInterval = Lambda(time>t1 & time<t2);

            SortedLambda_ = sort(DataLambdaInterval);
            SortedLambda = SortedLambda_(SortedLambda_>RatioMin&SortedLambda_<RatioMax);
            N = max(size(SortedLambda));

            x_exper = SortedLambda;
            y_exper_cdf = (1:1:N)/N;
            lambda = mean(SortedLambda);
            sigma = std(SortedLambda);
            x = RatioMin:0.005:RatioMax;
            y_cdf = cdf('Normal',x,lambda,sigma);
            y_pdf = pdf('Normal',x,lambda,sigma);
        end

    end


end