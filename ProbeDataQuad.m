classdef ProbeDataQuad < ProbeData & handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = protected, GetAccess = public)
        t;
        I1; 
        I2;
        U13;
        U43;
        CurrentNormilized;

        AreaRelativeI1 = 1;
        AreaRelativeI2 = 1;

        tgun;
        Igun;

        VAC_U13;
        VAC_U13_reduce;
        VAC_CurrentNormilized;
        VAC_CurrentNormilized_reduce;

        VAC_X_Fit;
        VAC_Y_Fit;

        VAC_IsCorrect;
        VAC_t1;
        VAC_t2;
        VAC_t;
        VAC_Isat;
        VAC_Isat_Dispersion;
        VAC_Er;
        VAC_Er_Dispersion;

        Temperature;
        Density;
        FloatPotential;

        Temperature_Dispersion;
        Density_Dispersion;
        FloatPotential_Dispersion;
        

        Goal;

    end


    methods (Access = public)

        function obj = ProbeDataQuad(ShotNumber,ProbeNumber,ElectronicBlockNumber,IndexPosition)
            
            obj = obj@ProbeData(ShotNumber,ProbeNumber,ElectronicBlockNumber,IndexPosition);
            
            if (nargin>0)   
                % Определяем состояния эксеримента:
                obj = WhatShotState(obj);
                % ProbeState = 0 or 1
                cd Shots\
                if obj.ProbeState > 0
                    I1_=importdata(['Ip1.',num2str(ShotNumber)]);
                    I2_=importdata(['Ipc1.',num2str(ShotNumber)]);
                    U13_=importdata(['Up1.',num2str(ShotNumber)]);
                    U43_=importdata(['Upd1.',num2str(ShotNumber)]);
                    obj.I1 = (I1_(1:end,2))/obj.AreaRelativeI1;
                    obj.I2 = (I2_(1:end,2))/obj.AreaRelativeI2;
                    obj.U13 = (U13_(1:end,2));
                    obj.U43 = (U43_(1:end,2));
                    obj.t = (I1_(1:end,1)/1000);
                    obj.CurrentNormilized = obj.I1./obj.I2;
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
                % Вычитание фона
                J1 = obj(s).I1;
                J2 = obj(s).I2;
                %V13 = obj(s).U13;
                V43 = obj(s).U43;
                Jgun = obj(s).Igun;
                I1_ground=sum(J1(1:50))/50;
                I2_ground=sum(J2(1:50))/50;
                U13_ground=0;
                U43_ground=sum(V43(1:50))/50;
                Igun_ground=sum(Jgun(1:50))/50;
                obj(s).I1=(obj(s).I1-I1_ground);
                obj(s).I2=(obj(s).I2-I2_ground);
                obj(s).U13=(obj(s).U13-U13_ground);
                obj(s).U43=(obj(s).U43-U43_ground);
                obj(s).Igun=(obj(s).Igun-Igun_ground);
                obj(s).CurrentNormilized = obj(s).I1./obj(s).I2;
                % Коррекция временных сигналов
                obj(s).t = obj(s).t - 0.1;


                % задание коэффициентов коррекции Амплитуд (эффекты блока электроники)
                % задание коэффициентов коррекции
                if obj(s).ElectronicBlockNumber == 1
                    r=12.5*10^-6;
                    r0=-0.5*10^-3;
                    k=6.5*10^-6;
                    k0=0.1*10^-3;
                elseif obj(s).ElectronicBlockNumber == 2
                end

                % Коррекция Амплитуд сигналов
                dI1=r*obj(s).U13;
                obj(s).I1 = obj(s).I1 - dI1+r0;
                dI2=k*obj(s).U13;
                obj(s).I2 = obj(s).I2 + dI2+k0;

                obj(s).IsSignalsCorrected = 1;
            end

            if (obj(s).PlasmaGunState == 1)&&(obj(s).IsSignalsCorrected == 0)
                % Вычитание фона
                Jgun = obj(s).Igun;
                Igun_ground=sum(Jgun(1:50))/50;
                obj(s).Igun=(obj(s).Igun-Igun_ground);
                % Коррекция временных сигналов
                obj(s).tgun = obj(s).tgun - 0.035;                
            end
            end
        end

        function SetAutoCalibration(obj)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
                if obj(s).IsSignalsCalibrated == 0
                    RatioInterval = [0 2];
                    CumMu = 1;
                    for j = 1:5
                        for i = 1:11
                            time(i) = 0.44+0.363*(i-1);
                            TimeInterval =[-0.01 0.01] + 0.44+0.363*(i-1);
                            [mu1(i), sigma1(i),t_interval,ratio_interval,x_Experimental,y_Experimental_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj(s),TimeInterval,RatioInterval);
                        end
                        mu(s) = mean(mu1(1:6));
                        sigma(s) = std(mu1(1:6));
                        obj(s).I1 = obj(s).I1/mu(s);
                        obj(s).I2 = obj(s).I2/1;
                        obj(s).CurrentNormilized = obj(s).I1./obj(s).I2;
                        CumMu = CumMu*mu(s);
                    end
                    obj(s).IsSignalsCalibrated = 1;
                    obj(s).AreaRelativeI1 = CumMu;
                    obj(s).AreaRelativeI2 = 1;
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
            goal = obj(s).Goal;
        end

        function DrawPlot(obj,limit_Y_1,limit_Y_2)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
            if obj(s).ProbeState == 1
                               
            f = figure(obj(s).ShotNumber);
            f.Color = [1 1 1];
            f.Position = [170 50 1200 730];

            tiledlayout(2,1)

            % Top plot
            ax1 = nexttile;
            p1 = plot(ax1,obj(s).t,obj(s).I1,obj(s).t,obj(s).I2);
            p1(1).Color  = [1 0 0];
            p1(2).Color  = [0 0 0];
            p1(1).LineWidth = 1;
            p1(2).LineWidth = 1;
            s1  = gca;
            s1.XLim = [-0.5 4.25];
            s1.YLim = limit_Y_1;
            % s1.YLim = [-0.1 0.3];
            % s1.XGrid = 'on';
            % s1.YGrid = 'on';
            s1.Title.String = 'Токи на зондовые проволоки 1 и 2';
            s1.XLabel.String = 't, мс';
            s1.YLabel.String = 'I, А';
            s1.FontSize = 16;
            s1.FontName = 'Times New Roman';
            s1.FontAngle = 'italic';
            % s1.Color = [1 1 0.5];

            % Bottom plot
            ax2 = nexttile;
            p2 = plot(ax2,obj(s).t,obj(s).U13,obj(s).t,obj(s).U43);
            p2(1).Color  = [0.00 0.45 0.74];
            p2(2).Color  = [0 0 0];
            p2(1).LineWidth = 1;
            p2(2).LineWidth = 1;

            s2  = gca;
            s2.XLim = [-0.5 4.25];
            s2.YLim = limit_Y_2;
            % s2.YLim = [-50 200];
            % s2.XGrid = 'on';
            % s2.YGrid = 'on';
            s2.Title.String = 'Напряжения между зондоыми проволоками 13 и 43';
            s2.XLabel.String = 't, мс';
            s2.YLabel.String = 'U, В';
            s2.FontSize = 16;
            s2.FontName = 'Times New Roman';
            s2.FontAngle = 'italic';
            % s2.Color = [1 1 0.5];
            else
                disp('There is no Shot in folder: D:\Евгений\The Class ProbeData\Shots');
            end
            end
        end

        function DrawCurrentNormilized(obj,limit_Y)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
                if obj(s).ProbeState == 1
                    f = figure(obj(s).ShotNumber);
                    f.Color = [1 1 1];
                    f.Position = [170 50 1200 730];
                    ax1 = nexttile;
                    p1 = plot(ax1,obj(s).t,obj(s).CurrentNormilized);
                    p1(1).Color  = [1 0 0];
                    %             p1(2).Color  = [0 0 0];
                    p1(1).LineWidth = 1;
                    %             p1(2).LineWidth = 1;
                    s1  = gca;
                    s1.XLim = [-0.5 4.25];
                    s1.YLim = limit_Y;
                    % s1.YLim = [-0.1 0.3];
                    % s1.XGrid = 'on';
                    % s1.YGrid = 'on';
                    s1.Title.String = 'Токи на зондовые проволоки 1 и 2';
                    s1.XLabel.String = 't, мс';
                    s1.YLabel.String = 'I, А';
                    s1.FontSize = 16;
                    s1.FontName = 'Times New Roman';
                    s1.FontAngle = 'italic';
                end
            end
        end

        function [mu, sigma, time] = GetCurrentNormilizedIonParthTimePoints(obj)
            RatioInterval = [0 2];
            for i =1:11
            time(i) = 0.44+0.363*(i-1);
            TimeInterval =[-0.01 0.01] + 0.44+0.363*(i-1);
            [mu(i), sigma(i),t_interval,ratio_interval,x_Experimental,y_Experimental_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj,TimeInterval,RatioInterval);
            end
        end

        function [mu, sigma] = GetCurrentNormilizedIonParth(obj)

            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
                RatioInterval = [0 2];
                for i =1:11
                    time(i) = 0.44+0.363*(i-1);
                    TimeInterval =[-0.01 0.01] + 0.44+0.363*(i-1);
                    [mu1(i), sigma1(i),t_interval,ratio_interval,x_Experimental,y_Experimental_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj(s),TimeInterval,RatioInterval);
                end
                mu(s) = mean(mu1(1:6));
                sigma(s) = std(mu1(1:6));
            end

        end

        function [Uexp,Jexp,t1,t2,tm,Er,Isat] = GetVac(obj)
            Npos = 23;
            T0pos = 2884-510*3;
            Tpos = 508;
            pos1 = T0pos-Tpos/2:Tpos/2:(T0pos-Tpos/2)+(Npos-1)*Tpos/2;
            pos2 = T0pos:Tpos/2:(T0pos)+(Npos-1)*Tpos/2;

            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
            if (obj(s).ProbeState == 1)&&(obj(s).IsSignalsCorrected == 1)&&(obj(s).IsSignalsCalibrated == 1)
            for j=1:Npos
            U_=obj(s).U13;I_=obj(s).CurrentNormilized;t_=obj(s).t;Isat_=obj(s).I2;
            Er_=obj(s).U43;
          
            Uexp(s,j) = {U_(pos1(j):pos2(j))};
            Jexp(s,j) = {I_(pos1(j):pos2(j))};
            Isat(s,j) = -mean(Isat_(pos1(j):pos2(j)));
            dIsat(s,j) = std(Isat_(pos1(j):pos2(j)));
            Er(s,j) = mean(Er_(pos1(j):pos2(j)));
            dEr(s,j) = std(Er_(pos1(j):pos2(j)));
            t1(j) = t_(pos1(j));t2(j) = t_(pos2(j)); tm(j) = (t1(j)+t2(j))/2;
            end
            end
            obj(s).VAC_U13 = Uexp(s,1:end);
            obj(s).VAC_CurrentNormilized = Jexp(s,1:end);
            obj(s).VAC_t1 = t1;
            obj(s).VAC_t2 = t2;
            obj(s).VAC_t = tm;
            obj(s).VAC_Isat = Isat(s,1:end);
            obj(s).VAC_Isat_Dispersion = dIsat(s,1:end);
            obj(s).VAC_Er = Er(s,1:end);
            obj(s).VAC_Er_Dispersion = dEr(s,1:end);
            end
        end

        function [Uexp_reduce,Jexp_reduce,X_Fit,Y_Fit,n,Te,U3,dn,dTe,dU3] = GetPlasmaParameters(obj,bounds_U,bounds_I)
            Npos = 23;
            ObjectNumber = max(size(obj));

            for s = 1:ObjectNumber
                if (obj(s).ProbeState == 1)&&(obj(s).IsSignalsCorrected == 1)&&(obj(s).IsSignalsCalibrated == 1)
                    cells_VAC_U13 = obj(s).VAC_U13;
                    cells_VAC_CurrentNormilized = obj(s).VAC_CurrentNormilized;
                    cells_VAC_Isat = obj(s).VAC_Isat;
                    cells_VAC_Isat_Dispersion = obj(s).VAC_Isat_Dispersion;

                    for j=1:Npos
                        Uexp = cells_VAC_U13{j};
                        Jexp = cells_VAC_CurrentNormilized{j};
                        J = cells_VAC_Isat(j);
                        Jstd = cells_VAC_Isat_Dispersion(j);
                        index_I = (bounds_U(1)<Uexp)&(Uexp<bounds_U(2));
                        index_U = (bounds_I(1)<Jexp)&(Jexp<bounds_I(2));

                        index = index_U&index_I;
                        Uexp_reduce(s,j) = {Uexp(index)};
                        Jexp_reduce(s,j) = {Jexp(index)};


                        X0 =  Uexp(index);
                        Y0 =  Jexp(index);
                        Te0=5;U20=-50;R=50;CirclesNumber=5;
                        [n_,Te_,U3_,dn_,dTe_,dU3_,etta_] = Plasma_Density_Circles(J,Jstd,Te0,U20,CirclesNumber,X0,Y0);
                        n(s,j)=n_;
                        Te(s,j)=Te_;
                        U3(s,j)=U3_;
                        dn(s,j)=dn_;
                        dTe(s,j)=dTe_;
                        dU3(s,j)=dU3_;
                        etta(s,j)=etta_;

                        X = (R/0.743)*sqrt(n_/Te_);
                        a1 = 0.3088*exp(-0.2225*X)+1.087;
                        b1 = -4.043*X.^0.03034+4.616;

                        u2 = -U3_+30;
                        x = -30:0.1:10;
                        y = ((x-u2+30)./(-u2)).^b1-42.8509*((u2/Te_).^-(b1)).*(1/a1).*exp((-u2+30)/Te_).*exp(x/Te_);
                        X_Fit(s,j) = {x};
                        Y_Fit(s,j) = {y};

                        
                    end

                else
                    Uexp_reduce(s,1:Npos) = nan;
                    Jexp_reduce(s,1:Npos) = nan;
                    X_Fit(s,1:Npos) = nan;
                    X_Fit(s,1:Npos) = nan;
                end

                obj(s).VAC_U13_reduce = Uexp_reduce(s,1:end);
                obj(s).VAC_CurrentNormilized_reduce = Jexp_reduce(s,1:end);

                obj(s).Density = n(s,1:end);
                obj(s).Temperature = Te(s,1:end);
                obj(s).FloatPotential = U3(s,1:end);
                obj(s).Density_Dispersion = dn(s,1:end);
                obj(s).Temperature_Dispersion = dTe(s,1:end);
                obj(s).FloatPotential_Dispersion = dU3(s,1:end);
                obj(s).VAC_X_Fit = X_Fit(s,1:end);
                obj(s).VAC_Y_Fit = Y_Fit(s,1:end);
            end

        end

        function DrawVAC(obj,X_lim,Y_lim,IntervalNumber)
            
            Npos = 23;
            cell_X = obj.VAC_U13_reduce;
            cell_Y = obj.VAC_CurrentNormilized_reduce;
            cell_X_Fit = obj.VAC_X_Fit;
            cell_Y_Fit = obj.VAC_Y_Fit;

            for j = 1:IntervalNumber(2)-IntervalNumber(1)+1

                f = figure(j);
                f.Color = [1 1 1];
                %f.Position = [170 50 1200 730];


                p = plot(cell_X{j},cell_Y{j},cell_X_Fit{j},cell_Y_Fit{j});
                p(1).LineStyle = 'none';
                p(1).Marker = '.';
                p(1).MarkerSize = 20;
                p(1).Color = 'k';
                p(2).Color = 'r';
                p(2).LineWidth = 2;

                s1  = gca;
                s1.XLim = X_lim;
                s1.YLim = Y_lim;
                %s1.XGrid = 'on';
                %s1.YGrid = 'on';
                t1 = round(obj.VAC_t1,2);
                t2 = round(obj.VAC_t2,2);
                title = ['NB',num2str(obj.ShotNumber),...
                                         '   t = ',num2str(t1(j)),'-',num2str(t2(j)),' мс'];
                x_label = 'U_{1-3}';
                y_label = 'I_{1} / I_{2}';
                s1.Title.String = title;
                s1.XLabel.String = x_label;
                s1.YLabel.String = y_label;
                s1.XLabel.FontWeight = 'bold';
                s1.YLabel.FontWeight = 'bold';
                s1.FontSize = 20;
                s1.FontName = 'Times New Roman';
                s1.FontAngle = 'italic';
                % s1.Color = [1 1 0.5];

                %                 leg{ObjectNumber} = 0;
                %                 for i = 1:ObjectNumber
                %                     q = ['NB',num2str(obj(i).timeInterval(1)),'-',num2str(obj(i).timeInterval(2)),...
                %                         '   t = ',num2str(obj(i).shotInterval(1)),'-',num2str(obj(i).shotInterval(2)),' мс'];
                %                     leg{i} = q;
                %                 end
                %
                %                 %legend(leg)
                %                 legend(leg{1})

            end

        end

        function [h,Er,dEr,t1,t2] = Get_VAC_Er(obj,IntervalNumber,VoltageShift)

            time1 = obj.VAC_t1; time2 = obj.VAC_t2;
            j=IntervalNumber;
            t1 = time1(j);t2 = time2(j);
            time = [t1 t2];
            ObjectNumber = max(size(obj));
            
            for s=1:ObjectNumber
                arr_Er = obj(s).VAC_Er;
                arr_dEr = obj(s).VAC_Er_Dispersion;
                arr_h = obj(s).ProbePositionAxisH;
                Er(s)=arr_Er(j)+VoltageShift;
                dEr(s)=arr_dEr(j);
                h(s)=arr_h;
            end

        end


    end



    methods (Access = private)

        function obj = WhatShotState(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
                        cd Shots\

            fileID_Ip1 = fopen(['Ip1.' num2str(obj.ShotNumber)]);
            fileID_Ipc1 = fopen(['Ipc1.' num2str(obj.ShotNumber)]);
            fileID_Up1 = fopen(['Up1.' num2str(obj.ShotNumber)]);
            fileID_Upd1 = fopen(['Upd1.' num2str(obj.ShotNumber)]);
            fileID_Igun = fopen(['Gun_Cath_Curr.',num2str(obj.ShotNumber)]);

            fileID_Probe = (fileID_Ip1 >= 0)&(fileID_Ipc1 >= 0)&(fileID_Up1 >= 0)&(fileID_Upd1 >= 0);
            fileID_Gun = (fileID_Igun >= 0);

            obj.ProbeState = fileID_Probe;
            obj.PlasmaGunState = fileID_Gun;
            obj.ShotState = (fileID_Gun&fileID_Probe);
            cd ../
        end

        function [mu, sigma,t_interval,ratio_interval,x_Experimental,y_Experimental_cdf, x, y_cdf, y_pdf] = CalculateCurrentRatio(obj,TimeInterval,RatioInterval)
            ObjectNumber = max(size(obj));
            time = obj(1).t;

            ratio = 0;
            for i=1:ObjectNumber
                ratio = ratio + obj(i).CurrentNormilized;
            end
            ratio=ratio/ObjectNumber;

            t1 = TimeInterval(1);           t2 = TimeInterval(2);           % Интервал времени
            RatioMin = RatioInterval(1);    RatioMax = RatioInterval(2);    % Интервал отнношений
            t_interval = time(time>t1 & time<t2);   ratio_interval = ratio(time>t1 & time<t2);


            SortedRatio_ = sort(ratio_interval);
            SortedRatio = SortedRatio_(SortedRatio_>RatioMin&SortedRatio_<RatioMax);  

            N = max(size(SortedRatio));
            RatioDistibution = (1:1:N)/N;

            x_Experimental = RatioDistibution;
            y_Experimental_cdf = SortedRatio;

            mu = mean(SortedRatio);
            sigma = std(SortedRatio);

            x = RatioMin:0.005:RatioMax;
            y_cdf = cdf('Normal',x,mu,sigma);
            y_pdf = pdf('Normal',x,mu,sigma);
        end

    end     
        

end


function J = Plasma_Current(n,Te,U2,R)
X = R./(0.743*sqrt(Te./n));
Y = -U2./Te;
a = 0.3088*exp(-0.2225*X)+1.087;
b = -4.043*X.^0.03034+4.616;

if (X>=4)&(X<=25)&(Y>=1)&(Y<=25)&(n>0)
I_plus_minus = a.*Y.^b;
else
I_plus_minus = nan;
end

e = 1.6*10^-19;
n0 = 10^20;%unit 10^14 cm^-3
mi = 1.675*10^-27;
S = 2*pi*(R*10^-6)*(10^-3);
Vs = sqrt(2*e*Te/mi);

    J = (1/sqrt(2*exp(1)))*e*n0.*n*S.*Vs.*sqrt(exp(1)/(2*pi)).*I_plus_minus;

end

function n = Plasma_Density(J,Te,U0,R,iteration)

% n - ед. 10^14 см^-3
% J - Амперы;
if (J>10^-4)
a0 = 1.087;
a1 = 0.3088;
a2 = 0.2225;
b0 = 4.616;
b1 = -4.043;
b2 = 0.03034;

e = 1.6*10^-19;
n0 = 10^20;%unit 10^14 cm^-3
mi = 1.675*10^-27;
S = 2*pi*(R*10^-6)*(10^-3);
Vs = sqrt(2*e*Te/mi);
k0=3.9198*10^-4;

%% Calculation
syms x
assume(x,'real')
assume(x>0)
assume(x<5)

eqn = (J == k0*sqrt(Te)*R*x*...
    (a1*exp(-a2*(R/0.743)*sqrt(x)/sqrt(Te))+a0)*...
    (-U0/Te)^(b1*(((R/0.743)^b2)*(sqrt(x)^b2))*(sqrt(Te)^-b2)+b0));

n = double(vpasolve(eqn, x));
X = R./(0.743*sqrt(Te./n));
Y = -U0./Te;

if iteration == 0
if (X>=4)&(X<=25)&(Y>=1)&(Y<=25)
else
n=nan;
end
end

else
n=nan;
end

end

function [n,Te,U3,dn,dTe,dU3,etta] = Plasma_Density_Circles(J,Jstd,Te,U2,CirclesNumber,U,I)

if (J>10^-4)&(J>Jstd)
% Calculate #1
R = 50;

n = Plasma_Density(J,Te,U2,R,1);
X = (R/0.743)*sqrt(n/Te);
a1 = 0.3088*exp(-0.2225*X)+1.087;
b1 = -4.043*X.^0.03034+4.616;

[fitresult, gof] = createFit_volt_ampere_characteristic(U, I, a1, b1);
eq = formula(fitresult);
coeffs = coeffvalues(fitresult);
ci = confint(fitresult,0.997);
for i=1:CirclesNumber-1
Te = coeffs(1);
U2 = -coeffs(2);
n = Plasma_Density(J,Te,U2,R,1);
X = (R/0.743)*sqrt(n/Te);
a1 = 0.3088*exp(-0.2225*X)+1.087;
b1 = -4.043*X.^0.03034+4.616;
[fitresult, gof] = createFit_volt_ampere_characteristic(U, I, a1, b1);

eq = formula(fitresult);
coeffs = coeffvalues(fitresult);
ci = confint(fitresult,0.997);
end


Te = coeffs(1);
U0 = coeffs(2);
n = Plasma_Density(J,Te,-U0,R,0);

if isnan(n)
    Te=nan;
    U3=nan;
    etta = nan;
    dn = nan;
    dTe = nan;
    dU3 = nan;
    A1 = nan;
    B1 = nan;

else
U3 = -U0+30;
etta = (U0-30)/Te;
dn = (Jstd/J)*n;
dTe = 0.5*(ci(2,1)-ci(1,1));
dU3 = 0.5*(ci(2,2)-ci(1,2));

X = (R/0.743)*sqrt(n/Te);
A1 = 0.3088*exp(-0.2225*X)+1.087;
B1 = -4.043*X.^0.03034+4.616;
end

else
n=nan;Te=nan;U3=nan;dn=nan;dTe=nan;dU3=nan;etta=nan;A1=nan;B1=nan;    
end

end

function [fitresult, gof] = createFit_volt_ampere_characteristic(U, I, a1, b1)
% b1=2.765*X^-0.052-2.182;
% a1=(0.49*exp(-0.306*X)+1.1)

[xData, yData] = prepareCurveData( U, I );

% Set up fittype and options.
% [Te,U0] - variable
% {a,b} - fixed
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1,30],...
               'Upper',[25,100],...
               'StartPoint',[5, 50]);
ft = fittype('((x-U0+30)/(-U0))^b-42.8509*((U0/Te)^-(b))*(1/a)*exp((-U0+30)/Te)*exp(x/Te)','problem',{'a','b'},'options',fo);
[fitresult, gof] = fit(xData, yData,ft,'problem',{a1,b1})

end

% Учет геометрических факторов зондов (коррекция токов)
% R0 = 50; L0 = 1000; L1 = 500; R1 = 100;
% S0 = 2*pi*R0*L0+pi*R0^2;
% S1 = 2*pi*R0*L1+(2*pi*R1^2)*(1+sqrt(1-(R0/R1)^2));
% S10 = S1/S0;
% S10 = 1;