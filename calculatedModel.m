classdef calculatedModel < handle
    
    %{
    1/5/21 - Chris Mills
    To be used with the Mills Senior Thesis; ensembleAnalysis3. Purely for 
    convenience of analyzing a number of different models.
    %}
    
    properties (Access = public) %public for easy editing
        %%%%%% Model properties %%%%%
        depths;     %array of model-space parameter values
        rhos;       %array of model-space parameter values
        %%%%%% Data properties %%%%%%%
        y;          %array of data-space output values
        misfit;     %real number
        %%%%%% Plotting properties %%%%%
        color;      %RGB triplet or char or hexadecimal or string
        lineStyle;  %char using standard matlab conventions
        displayName;      %Title/name of the model
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    methods
        
%%%%%%%%%%%%%%%% Core Functions %%%%%%%%%%%%%%%%%%%%%%

        %Constructor%
        function obj = calculatedModel(inDepths,inRhos,inY,dataY,...
                colorChoice,lineStyle,title)
            %inDepths, inRhos are arrays containing depths and 
            %resistivities associated with this model; lineStyle,
            %colorChoice,and title are for graphing and follow matlab 
            %plotting conventions for LineStyle, Color, and DisplayName 
            %respectively
            obj.depths = inDepths;
            obj.rhos = inRhos;
            obj.setY(inY,dataY);
            obj.color = colorChoice;
            obj.displayName = title;
            obj.lineStyle = lineStyle;
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
    end
end