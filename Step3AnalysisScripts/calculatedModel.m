classdef calculatedModel < handle
    
    %{
    1/5/21 - Chris Mills
    To be used with the Mills Senior Thesis; ensembleAnalysisA. Purely for 
    convenience of analyzing a number of different models.
    %}
    
    properties (Access = public) %public for easy editing
        %%%%%% Model properties %%%%%
        depths;     %array of model-space parameter values
        rhos;       %array of model-space parameter values
        %%%%%% Data properties %%%%%%%
        y;          %array of data-space output values
        misfit;     %real number: data-space misfit
        residual;   %
        wre2n;      %weighted relative error in the 2 norm
        %%%%%% Plotting properties %%%%%
        color;      %RGB triplet or char or hexadecimal or string
        lineStyle;  %char using standard matlab conventions
        displayName;%Title/name of the model as it'll appear in fig legends
        lineWidth;  %Defaults to 1.0, change it after if you want
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    methods
        
%%%%%%%%%%%%%%%% Core Functions %%%%%%%%%%%%%%%%%%%%%%

        %Constructor%
        function obj = calculatedModel(inDepths,inRhos,inY,dataY,...
                colorChoice,lineStyle,title)
            %inDepths, inRhos are arrays containing depths and 
            %resistivities associated with this model. lineStyle,
            %colorChoice,and title are for graphing and follow matlab 
            %plotting conventions for LineStyle, Color, and DisplayName 
            obj.depths = inDepths;
            obj.rhos = inRhos;
            obj.setY(inY,dataY);
            obj.color = colorChoice;
            obj.displayName = title;
            obj.lineStyle = lineStyle;
            obj.lineWidth = 1.0;
        end
        
%%%%%%%%%%%%%%% Set functions %%%%%%%%%%%%%%%%%%%%%

        function setY(obj,inY,dataY)
            %inY is the y-values in data space associated with this model 
            %(in our case, resistivities as a fxn of electrode spacings;
            %dataY is the 'actual' data from field measurements (or
            %synth-generated), used to calculate misfit. Both are arrays of
            %the same size
            obj.y = inY; 
            obj.setMisfit(dataY);
        end
        
        function setMisfit(obj,dataY)
            %see setY fxn
            obj.misfit = norm(obj.y - dataY);
        end
        
        function setWRE2N(obj,data)
            obj.residual = obj.y - data.y;
            obj.wre2n = sqrt(obj.residual'*inv(data.Cd)*obj.residual)/...
                sqrt(data.y'*inv(data.Cd)*data.y);
        end
    end
end