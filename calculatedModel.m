classdef calculatedModel < handle
    
    properties (Access = public)
        %%%%%% Model properties %%%%%
        depths;     %array
        rhos;       %array
        %%%%%% Data properties %%%%%%%
        y;          %array
        misfit;     %real number
        %%%%%% Plotting properties %%%%%
        color;      %RGB triplet (can be input as char)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    methods
        
%%%%%%%%%%%%%%%% Core Functions %%%%%%%%%%%%%%%%%%%%%%

        %Constructor%
        function obj = calculatedModel(inDepths,inRhos,...
                colorChoice)
            obj.depths = inDepths;
            obj.rhos = inRhos;
            obj.y = 0;
            obj.misfit = 0;
            obj.setColor(colorChoice);
        end
        
%%%%%%%%%%%%%%% Set functions %%%%%%%%%%%%%%%%%%%%%

        function setY(obj,inY,dataY)
            obj.y = inY;
            obj.setMisfit(dataY);
        end
        
        function setMisfit(obj,dataY)
            obj.misfit = norm(obj.y - dataY);
        end
        
        function setColor(obj,colorChoice)
            if isstring(colorChoice)
                switch colorChoice
                    case 'y'
                        obj.color = [1 1 0];
                    case 'm'
                        obj.color = [1 0 1];
                    case 'c'
                        obj.color = [0 1 1];
                    case 'r'
                        obj.color = [1 0 0];
                    case 'g'
                        obj.color = [0 1 0];
                    case 'b'
                        obj.color = [0 0 1];
                    case 'w'
                        obj.color = [1 1 1];
                    case 'k'
                        obj.color = [0 0 0];
                end
            else
                obj.color = colorChoice;
            end
        end       
    end
end