classdef ProbeDataDistribution < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        h;
        Data;
        DataDispersion;
        h1;
        h2;
        Data1;
        Data2;
        DataDispersion1;
        DataDispersion2;
        H_Min;

        h12;
        R12;
        Data12;
        DataDispersion12;

        timeInterval;
        shotInterval;

    end

    methods (Access = public)

        function obj = ProbeDataDistribution(h,Data,DataDispersion,timeInterval,shotInterval)
            obj.h = h;
            obj.Data = Data;
            obj.DataDispersion = DataDispersion;
            obj.timeInterval = timeInterval;
            obj.shotInterval = shotInterval;
        end

        function GetSymmetricDistribution(obj,h_DeliteRange)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
            [h_Union, Data_Union, DataDispersion_Union] = UnionData(obj(s).h,obj(s).Data,obj(s).DataDispersion);
            [h_AfterDelite,Data_AfterDelite,DataDispersion_AfterDelite] = DeliteData(h_Union, Data_Union,DataDispersion_Union,h_DeliteRange);
            [H_Min_] = Find_H_MirrorPositionSymmetric(h_AfterDelite,Data_AfterDelite,DataDispersion_AfterDelite);

            h1_ = h_Union; Data1_ = Data_Union;DataDispersion1_ = DataDispersion_Union;
            [h2_,Data2_,DataDispersion2_] = MirrorDataSymmetic(h1_,Data1_,DataDispersion1_,H_Min_);
            obj(s).h1 = h1_;
            obj(s).h2 = h2_;
            obj(s).h12 = [h1_,h2_];
            obj(s).R12 = [h1_,h2_]-H_Min_;
            obj(s).Data1 = Data1_;
            obj(s).Data2 = Data2_;
            obj(s).Data12 = [Data1_,Data2_];
            obj(s).DataDispersion1 = DataDispersion1_;
            obj(s).DataDispersion2 = DataDispersion2_;
            obj(s).DataDispersion12 = [DataDispersion1_,DataDispersion2_];
            obj(s).H_Min = H_Min_;

            end
        end

        function GetAntiSymmetricDistribution(obj,h_DeliteRange)
            ObjectNumber = max(size(obj));
            for s = 1:ObjectNumber
            [h_Union, Data_Union, DataDispersion_Union] = UnionData(obj(s).h,obj(s).Data,obj(s).DataDispersion);
            [h_AfterDelite,Data_AfterDelite,DataDispersion_AfterDelite] = DeliteData(h_Union, Data_Union,DataDispersion_Union,h_DeliteRange);
            [H_Min_] = Find_H_MirrorPositionAntiSymmetric(h_AfterDelite,Data_AfterDelite,DataDispersion_AfterDelite);

            h1_ = h_Union; Data1_ = Data_Union;DataDispersion1_ = DataDispersion_Union;
            [h2_,Data2_,DataDispersion2_] = MirrorDataAntiSymmetic(h1_,Data1_,DataDispersion1_,H_Min_);
            obj(s).h1 = h1_;
            obj(s).h2 = h2_;
            obj(s).Data1 = Data1_;
            obj(s).Data2 = Data2_;
            obj(s).Data12 = [Data1_,Data2_];
            obj(s).R12 = [h1_,h2_]-H_Min_;
            obj(s).DataDispersion1 = DataDispersion1_;
            obj(s).DataDispersion2 = DataDispersion2_;
            obj(s).DataDispersion12 = [DataDispersion1_,DataDispersion2_];
            obj(s).H_Min = H_Min_;
            end
        end
  
        function DrawDistribution(obj,X_lim,Y_lim,FigureNumber,title,x_label,y_label,plotType)

            colors{1} = [0 0 0];
            colors{2} = [1 0 0];
            colors{3} = [0.00,0.45,0.74];
            colors{4} = [0.47,0.67,0.19];

            ObjectNumber = max(size(obj));

            f = figure(FigureNumber);
            f.Color = [1 1 1];
%           f.Position = [170 50 1200 730];
            for s = 1:ObjectNumber
                if plotType == 0
                plot(obj(s).h, obj(s).Data,'LineStyle','none','Marker','.','MarkerSize',25,'Color',colors{s});
                elseif plotType == 1
                semilogy(obj(s).h, obj(s).Data,'LineStyle','none','Marker','.','MarkerSize',25,'Color',colors{s});
                end
                hold on;
            end
            hold off;

            s1  = gca;
            s1.XLim = X_lim;
            s1.YLim = Y_lim;
            %s1.XGrid = 'on';
            %s1.YGrid = 'on';
            s1.Title.String = title;
            s1.XLabel.String = x_label;
            s1.YLabel.String = y_label;
            s1.XLabel.FontWeight = 'bold';
            s1.YLabel.FontWeight = 'bold';
            s1.FontSize = 20;
            s1.FontName = 'Times New Roman';
            s1.FontAngle = 'italic';
            % s1.Color = [1 1 0.5];

            leg{ObjectNumber} = 0;
            for i = 1:ObjectNumber
                q = ['NB',num2str(obj(i).timeInterval(1)),'-',num2str(obj(i).timeInterval(2)),...
                    '   t = ',num2str(obj(i).shotInterval(1)),'-',num2str(obj(i).shotInterval(2)),' мс'];
                leg{i} = q;
            end

            legend(leg)

        end

        function DrawDistributionFit(obj,X_lim,Y_lim,SmoothingParam,title,x_label,y_label,figNumber)

            colors{1} = [0 0 0];
            colors{2} = [1 0 0];
            colors{3} = [0.00,0.45,0.74];
            colors{4} = [0.47,0.67,0.19];

            ObjectNumber = max(size(obj));

            for s = 1:ObjectNumber
                [fitresult, gof, xData, yData ] = createFit(obj(s).R12, obj(s).Data12,SmoothingParam);

                f = figure(s+figNumber-1);
                f.Color = [1 1 1];
                %           f.Position = [170 50 1200 730];

                p = plot( fitresult, xData, yData );
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
                s1.Title.String = title;
                s1.XLabel.String = x_label;
                s1.YLabel.String = y_label;
                s1.XLabel.FontWeight = 'bold';
                s1.YLabel.FontWeight = 'bold';
                s1.FontSize = 20;
                s1.FontName = 'Times New Roman';
                s1.FontAngle = 'italic';
                % s1.Color = [1 1 0.5];

                leg{ObjectNumber} = 0;

                for i = 1:ObjectNumber
                    t1 = round(obj(i).timeInterval(1),2);
                    t2 = round(obj(i).timeInterval(2),2);
                    q = ['NB',num2str(obj(i).shotInterval(1)),'-',num2str(obj(i).shotInterval(2)),...
                        '   t = ',num2str(t1),'-',num2str(t2),' мс'];
                    leg{i} = q;
                end

                %           legend(leg)
                legend(leg{1})
            

            end
           
        end

        function DrawErrorBarDistribution(obj,X_lim,Y_lim,FigureNumber,title,x_label,y_label)

            colors{1} = [0 0 0];
            colors{2} = [1 0 0];
            colors{3} = [0.00,0.45,0.74];
            colors{4} = [0.47,0.67,0.19];

            ObjectNumber = max(size(obj));

            f = figure(FigureNumber);
            f.Color = [1 1 1];
%           f.Position = [170 50 1200 730];
            for s = 1:ObjectNumber
                errorbar(obj(s).h, obj(s).Data, obj(s).DataDispersion,'LineStyle','none','Marker','.','MarkerSize',25,'Color',colors{s});
                hold on;
            end
            hold off;

            s1  = gca;
            s1.XLim = X_lim;
            s1.YLim = Y_lim;
            %s1.XGrid = 'on';
            %s1.YGrid = 'on';
            s1.Title.String = title;
            s1.XLabel.String = x_label;
            s1.YLabel.String = y_label;
            s1.XLabel.FontWeight = 'bold';
            s1.YLabel.FontWeight = 'bold';
            s1.FontSize = 20;
            s1.FontName = 'Times New Roman';
            s1.FontAngle = 'italic';
            % s1.Color = [1 1 0.5];

            leg{ObjectNumber} = 0;
            for i = 1:ObjectNumber
                q = ['NB',num2str(obj(i).timeInterval(1)),'-',num2str(obj(i).timeInterval(2)),...
                    '   t = ',num2str(obj(i).shotInterval(1)),'-',num2str(obj(i).shotInterval(2)),' мс'];
                leg{i} = q;
            end

            legend(leg)

        end

    end


end


function [h_Union, Data_Union, DataDispersion_Union,NumberUnionData] = UnionData(h,Data,DataDispersion)
h_Union = unique(h);
column = max(size(h));
column_Union = max(size(h_Union));
NumberUnionData = column - column_Union;
Data_Union(column_Union) = 0;
DataDispersion_Union(column_Union) = 0;
for i = 1:column_Union
    index = (h==h_Union(i));
    Data_Union(i)=mean(Data(index));
    DataDispersion_Union(i)=mean(DataDispersion(index));
end
end

function [h_AfterDelite,Data_AfterDelite,DataDispersion_AfterDelite] = DeliteData(h,Data,DataDispersion,h_DeliteRange)
index = (h<=h_DeliteRange(1))|(h>=h_DeliteRange(2));
h_AfterDelite = h(index);
Data_AfterDelite = Data(index);
DataDispersion_AfterDelite = DataDispersion(index);
end

function [h_Mirror,Data_Mirror,DataDispersion_Mirror] = MirrorDataSymmetic(h,Data,DataDispersion,H_MirrorPosition)
h_Mirror = 2*H_MirrorPosition-h;
Data_Mirror = Data;
DataDispersion_Mirror = DataDispersion;
end

function [h_Mirror,Data_Mirror,DataDispersion_Mirror] = MirrorDataAntiSymmetic(h,Data,DataDispersion,H_MirrorPosition)
h_Mirror = 2*H_MirrorPosition-h;
Data_Mirror = -Data;
DataDispersion_Mirror = DataDispersion;
end

function [H_Min,Chi_Square_Min,H_MirrorPosition_Variate,Chi_Square] = Find_H_MirrorPositionSymmetric(h,Data,DataDispersion)
hmin = min(h);
hmax = max(h);
step = 0.1;
H_MirrorPosition_Variate = hmin:step:hmax;
n = max(size(H_MirrorPosition_Variate));

Chi_Square(n) = 0;
for i = 1:n
    [h_Mirror,Data_Mirror] = MirrorDataSymmetic(h,Data,DataDispersion,H_MirrorPosition_Variate(i));
    Data_Extrapolation = interp1(h,Data,h_Mirror,'linear','extrap');
    Chi_Square(i) = sum((Data_Extrapolation-Data_Mirror).^2);
end
% Chi_Square = smooth(Chi_Square,100);
[value,index]=min(Chi_Square);
Chi_Square_Min = value;
H_Min = H_MirrorPosition_Variate(index);
end

function [H_Min,Chi_Square_Min,H_MirrorPosition_Variate,Chi_Square] = Find_H_MirrorPositionAntiSymmetric(h,Data,DataDispersion)
hmin = min(h);
hmax = max(h);
step = 0.1;
H_MirrorPosition_Variate = hmin:step:hmax;
n = max(size(H_MirrorPosition_Variate));

Chi_Square(n) = 0;
for i = 1:n
    [h_Mirror,Data_Mirror] = MirrorDataAntiSymmetic(h,Data,DataDispersion,H_MirrorPosition_Variate(i));
    Data_Extrapolation = interp1(h,Data,h_Mirror,'linear','extrap');
    Chi_Square(i) = sum((Data_Extrapolation-Data_Mirror).^2);
end
% Chi_Square = smooth(Chi_Square,100);
[value,index]=min(Chi_Square);
Chi_Square_Min = value;
H_Min = H_MirrorPosition_Variate(index);
end

function [fitresult, gof,  xData, yData] = createFit(X, Y, SmoothingParam)
%CREATEFIT3(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : X
%      Y Output: Y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 23-Aug-2022 12:41:06
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( X, Y);

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = SmoothingParam;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end