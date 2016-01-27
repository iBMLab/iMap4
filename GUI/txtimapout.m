function txtimapout(MeasureM,conditiontmp,CondName,opt)
% text file output
%
if opt==1 % overall descriptive result
    if isunix
        fileID = fopen(['Descriptive_STAT/Info_all.txt'],'w');
    else
        fileID = fopen(['Descriptive_STAT\Info_all.txt'],'w');
    end
    fprintf(fileID,'%s \n','Thank you for using iMap4. This .txt file contains the descriptive');
    fprintf(fileID,'%s \n','statistic of your data input.');
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','<Mean_Intensity_Map.png> is the mean fixation map of the whole dataset.');
    fprintf(fileID,'%s \n','<Eye_Mvt_Measure_Dist.png> is the histogram and the fitted distribution of the following eye');
    fprintf(fileID,'%s \n','movement measurement; these measurements are calculated per trial or');
    fprintf(fileID,'%s \n','condition depending on the chosen option:');
elseif opt==2 % condition or joint condition
    if isunix
        fileID = fopen(['Descriptive_STAT/',CondName,'/','Info_' CondName '.txt'],'w');
    else
        fileID = fopen(['Descriptive_STAT\',CondName,'\','Info_' CondName '.txt'],'w');
    end
    fprintf(fileID,'%s \n','Thank you for using iMap4. This .txt file contains the descriptive');
    fprintf(fileID,'%s \n','statistic of your data input.');
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n',['<', CondName, '_Mean_Intensity_Map.png> is the mean fixation map of each level in condition <' CondName, '>.']);
    fprintf(fileID,'%s \n',['<', CondName, '_Boxplot.png> is the boxplot of the following eye movement measurement for each']);
    fprintf(fileID,'%s \n',['level in condition <' CondName '>; these measurements are calculated per']);
    fprintf(fileID,'%s \n','trial or condition depending on the chosen option:');
end
fprintf(fileID,'%s \n','     ');
fprintf(fileID,'%s \n',' - Number of fixation');
fprintf(fileID,'%s \n',' - Sum of fixation duration (Total viewing time)');
fprintf(fileID,'%s \n',' - Mean fixation duration');
fprintf(fileID,'%s \n',' - Total path length (total eye movement path length in pixel)');
fprintf(fileID,'%s \n',' - Mean path length');
fprintf(fileID,'%s \n','     ');
fprintf(fileID,'%s \n','     ');
fprintf(fileID,'%s \n','Importantly, iMap4 fits a gamma distribution on all above measurement.');
fprintf(fileID,'%s \n','We believe this is appropriate and more informative than a Gaussian');
fprintf(fileID,'%s \n','assumption given that eye fixation should be considered as a Poisson');
fprintf(fileID,'%s \n','process. Thus, its parameters could be fit with a Gamma distribution.');
fprintf(fileID,'%s \n','Fixation duration and path length could both be considered some how as');
fprintf(fileID,'%s \n','the waiting time of the Poisson process, thus should also be appropriate');
fprintf(fileID,'%s \n','to fit with a gamma distribution.');

FixNum=MeasureM.FixNum;
sumFixDur=MeasureM.sumFixDur;
meanFixDur=MeasureM.meanFixDur;
totalPathLength=MeasureM.totalPathLength;
meanPathLength=MeasureM.meanPathLength;
if opt==2 % condition or joint condition
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n',['The following part is the number output of each level in condition <',CondName,'>.']);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','     ');
    
    uniquecondi=unique(conditiontmp);
    lengthcondi=length(uniquecondi);
    for icc=1:lengthcondi
        levelName=uniquecondi(icc);
        try
            idxtmp=conditiontmp==uniquecondi(icc);
        catch
            idxtmp=strcmp(conditiontmp,uniquecondi(icc));
        end
        try
            pd1 = fitdist(FixNum(idxtmp),'gamma');
            alpha1=pd1.a;beta1=pd1.b;paramci1=pd1.paramci;
        catch
            alpha1=NaN;beta1=NaN;paramci1=[NaN,NaN];
        end
        
        try
            pd2 = fitdist(sumFixDur(idxtmp),'gamma');
            alpha2=pd2.a;beta2=pd2.b;paramci2=pd2.paramci;
        catch
            alpha2=NaN;beta2=NaN;paramci2=[NaN,NaN];
        end
        
        try
            pd3 = fitdist(meanFixDur(idxtmp),'gamma');
            alpha3=pd3.a;beta3=pd3.b;paramci3=pd3.paramci;
        catch
            alpha3=NaN;beta3=NaN;paramci3=[NaN,NaN];
        end
        
        try
            pd4 = fitdist(totalPathLength(idxtmp),'gamma');
            alpha4=pd4.a;beta4=pd4.b;paramci4=pd4.paramci;
        catch
            alpha4=NaN;beta4=NaN;paramci4=[NaN,NaN];
        end
        
        try
            pd5 = fitdist(meanPathLength(idxtmp),'gamma');
            alpha5=pd5.a;beta5=pd5.b;paramci5=pd5.paramci;
        catch
            alpha5=NaN;beta5=NaN;paramci5=[NaN,NaN];
        end
        
        fprintf(fileID,'%s \n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fileID,'%s \n',['   Level Name: <' char(levelName) '>']);
        fprintf(fileID,'%s \n',['   Total number of observation ' num2str(sum(idxtmp==1))]);
        fprintf(fileID,'%s \n','     ');
        fprintf(fileID,'%s \n','    Number of fixation:');
        fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
        fprintf(fileID,'%s \n',['       ',num2str(quantile(FixNum(idxtmp),[0 0.25 0.50 0.75 1]))]);
        fprintf(fileID,'%s \n',['       Mean <' num2str(mean(FixNum(idxtmp))) '>; Std. Deviation <' num2str(std(FixNum(idxtmp))) '> ']);
        fprintf(fileID,'%s \n','      under Gamma distribution:');
        fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd1)) '>; Std. Deviation <' num2str(std(pd1)) '>']);
        fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha1) '>, [' num2str(paramci1(:,1)') '];']);
        fprintf(fileID,'%s \n',['        beta = <' num2str(beta1) '>, [' num2str(paramci1(:,2)') '];']);
        fprintf(fileID,'%s \n','     ');
        fprintf(fileID,'%s \n','    Sum of fixation duration:');
        fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
        fprintf(fileID,'%s \n',['       ',num2str(quantile(sumFixDur(idxtmp),[0 0.25 0.50 0.75 1]))]);
        fprintf(fileID,'%s \n',['       Mean <' num2str(mean(sumFixDur(idxtmp))) '>; Std. Deviation <' num2str(std(sumFixDur(idxtmp))) '> ']);
        fprintf(fileID,'%s \n','      under Gamma distribution:');
        fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd2)) '>; Std. Deviation <' num2str(std(pd2)) '>']);
        fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha2) '>, [' num2str(paramci2(:,1)') '];']);
        fprintf(fileID,'%s \n',['        beta = <' num2str(beta2) '>, [' num2str(paramci2(:,2)') '];']);
        fprintf(fileID,'%s \n','     ');
        fprintf(fileID,'%s \n','    Mean fixation duration:');
        fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
        fprintf(fileID,'%s \n',['       ',num2str(quantile(meanFixDur(idxtmp),[0 0.25 0.50 0.75 1]))]);
        fprintf(fileID,'%s \n',['       Mean <' num2str(mean(meanFixDur(idxtmp))) '>; Std. Deviation <' num2str(std(meanFixDur(idxtmp))) '> ']);
        fprintf(fileID,'%s \n','      under Gamma distribution:');
        fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd3)) '>; Std. Deviation <' num2str(std(pd3)) '>']);
        fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha3) '>, [' num2str(paramci3(:,1)') '];']);
        fprintf(fileID,'%s \n',['        beta = <' num2str(beta3) '>, [' num2str(paramci3(:,2)') '];']);
        fprintf(fileID,'%s \n','     ');
        fprintf(fileID,'%s \n','    Total path length:');
        fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
        fprintf(fileID,'%s \n',['       ',num2str(quantile(totalPathLength(idxtmp),[0 0.25 0.50 0.75 1]))]);
        fprintf(fileID,'%s \n',['       Mean <' num2str(mean(totalPathLength(idxtmp))) '>; Std. Deviation <' num2str(std(totalPathLength(idxtmp))) '> ']);
        fprintf(fileID,'%s \n','      under Gamma distribution:');
        fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd4)) '>; Std. Deviation <' num2str(std(pd4)) '>']);
        fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha4) '>, [' num2str(paramci4(:,1)') '];']);
        fprintf(fileID,'%s \n',['        beta = <' num2str(beta4) '>, [' num2str(paramci4(:,2)') '];']);
        fprintf(fileID,'%s \n','     ');
        fprintf(fileID,'%s \n','    Mean path length:');
        fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
        fprintf(fileID,'%s \n',['       ',num2str(quantile(meanPathLength(idxtmp),[0 0.25 0.50 0.75 1]))]);
        fprintf(fileID,'%s \n',['       Mean <' num2str(mean(meanPathLength(idxtmp))) '>; Std. Deviation <' num2str(std(meanPathLength(idxtmp))) '> ']);
        fprintf(fileID,'%s \n','      under Gamma distribution:');
        fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd5)) '>; Std. Deviation <' num2str(std(pd5)) '>']);
        fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha5) '>, [' num2str(paramci5(:,1)') '];']);
        fprintf(fileID,'%s \n',['        beta = <' num2str(beta5) '>, [' num2str(paramci5(:,2)') '];']);
        fprintf(fileID,'%s \n','     ');
        
    end
else
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','The following part is the number output of the whole experiment');
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','     ');
    
    try
        pd1 = fitdist(FixNum,'gamma');
        alpha1=pd1.a;beta1=pd1.b;paramci1=pd1.paramci;
    catch
        alpha1=NaN;beta1=NaN;paramci1=[NaN,NaN];
    end
    
    try
        pd2 = fitdist(sumFixDur,'gamma');
        alpha2=pd2.a;beta2=pd2.b;paramci2=pd2.paramci;
    catch
        alpha2=NaN;beta2=NaN;paramci2=[NaN,NaN];
    end
    
    try
        pd3 = fitdist(meanFixDur,'gamma');
        alpha3=pd3.a;beta3=pd3.b;paramci3=pd3.paramci;
    catch
        alpha3=NaN;beta3=NaN;paramci3=[NaN,NaN];
    end
    
    try
        pd4 = fitdist(totalPathLength,'gamma');
        alpha4=pd4.a;beta4=pd4.b;paramci4=pd4.paramci;
    catch
        alpha4=NaN;beta4=NaN;paramci4=[NaN,NaN];
    end
    
    try
        pd5 = fitdist(meanPathLength,'gamma');
        alpha5=pd5.a;beta5=pd5.b;paramci5=pd5.paramci;
    catch
        alpha5=NaN;beta5=NaN;paramci5=[NaN,NaN];
    end
    
    fprintf(fileID,'%s \n',['   Total number of observation ' num2str(length(FixNum))]);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','    Number of fixation:');
    fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
    fprintf(fileID,'%s \n',['       ',num2str(quantile(FixNum,[0 0.25 0.50 0.75 1]))]);
    fprintf(fileID,'%s \n',['       Mean <' num2str(mean(FixNum)) '>; Std. Deviation <' num2str(std(FixNum)) '> ']);
    fprintf(fileID,'%s \n','      under Gamma distribution:');
    fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd1)) '>; Std. Deviation <' num2str(std(pd1)) '>']);
    fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha1) '>, [' num2str(paramci1(:,1)') '];']);
    fprintf(fileID,'%s \n',['        beta = <' num2str(beta1) '>, [' num2str(paramci1(:,2)') '];']);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','    Sum of fixation duration:');
    fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
    fprintf(fileID,'%s \n',['       ',num2str(quantile(sumFixDur,[0 0.25 0.50 0.75 1]))]);
    fprintf(fileID,'%s \n',['       Mean <' num2str(mean(sumFixDur)) '>; Std. Deviation <' num2str(std(sumFixDur)) '> ']);
    fprintf(fileID,'%s \n','      under Gamma distribution:');
    fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd2)) '>; Std. Deviation <' num2str(std(pd2)) '>']);
    fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha2) '>, [' num2str(paramci2(:,1)') '];']);
    fprintf(fileID,'%s \n',['        beta = <' num2str(beta2) '>, [' num2str(paramci2(:,2)') '];']);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','    Mean fixation duration:');
    fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
    fprintf(fileID,'%s \n',['       ',num2str(quantile(meanFixDur,[0 0.25 0.50 0.75 1]))]);
    fprintf(fileID,'%s \n',['       Mean <' num2str(mean(meanFixDur)) '>; Std. Deviation <' num2str(std(meanFixDur)) '> ']);
    fprintf(fileID,'%s \n','      under Gamma distribution:');
    fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd3)) '>; Std. Deviation <' num2str(std(pd3)) '>']);
    fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha3) '>, [' num2str(paramci3(:,1)') '];']);
    fprintf(fileID,'%s \n',['        beta = <' num2str(beta3) '>, [' num2str(paramci3(:,2)') '];']);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','    Total path length:');
    fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
    fprintf(fileID,'%s \n',['       ',num2str(quantile(totalPathLength,[0 0.25 0.50 0.75 1]))]);
    fprintf(fileID,'%s \n',['       Mean <' num2str(mean(totalPathLength)) '>; Std. Deviation <' num2str(std(totalPathLength)) '> ']);
    fprintf(fileID,'%s \n','      under Gamma distribution:');
    fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd4)) '>; Std. Deviation <' num2str(std(pd4)) '>']);
    fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha4) '>, [' num2str(paramci4(:,1)') '];']);
    fprintf(fileID,'%s \n',['        beta = <' num2str(beta4) '>, [' num2str(paramci4(:,2)') '];']);
    fprintf(fileID,'%s \n','     ');
    fprintf(fileID,'%s \n','    Mean path length:');
    fprintf(fileID,'%s \n','        min     1st quartile    median     3rd quartile     max');
    fprintf(fileID,'%s \n',['       ',num2str(quantile(meanPathLength,[0 0.25 0.50 0.75 1]))]);
    fprintf(fileID,'%s \n',['       Mean <' num2str(mean(meanPathLength)) '>; Std. Deviation <' num2str(std(meanPathLength)) '> ']);
    fprintf(fileID,'%s \n','      under Gamma distribution:');
    fprintf(fileID,'%s \n',['        Mean <' num2str(mean(pd5)) '>; Std. Deviation <' num2str(std(pd5)) '>']);
    fprintf(fileID,'%s \n',['        alpha = <' num2str(alpha5) '>, [' num2str(paramci5(:,1)') '];']);
    fprintf(fileID,'%s \n',['        beta = <' num2str(beta5) '>, [' num2str(paramci5(:,2)') '];']);
end

fclose(fileID);
