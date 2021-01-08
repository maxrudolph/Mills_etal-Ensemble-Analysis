function doSaving(filename,results,data,measure,options,forwardModel,pBounds)
%Matlab won't allow you to put a save statement in a parfor loop, but if 
%you put it in a function then its okay for some reason. 
     save(filename,'results','data','measure','options','forwardModel',...
         'pBounds','-v7.3'); %-v7.3 allows for saving of large files
end
