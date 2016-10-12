function imapLMMreport(StatMap_c,varargin)
if nargin>1
    % a 4d matrix with (Npredictor,beta[95%CI],xSize,ySize)
    betamaps  = varargin{1};
    plotcondi = 1;
else
    plotcondi = 0;
end
%%
clc
label = StatMap_c.label;
Pmask = StatMap_c.Pmask;
Fmap  = StatMap_c.map;
dof   = StatMap_c.df;
if isfield(StatMap_c,'beta')
    betaall   = StatMap_c.beta;
    betaCIall = StatMap_c.betaCI;
end
for ilabel = 1:length(label)
    Pmasktmp = squeeze(Pmask(ilabel,:,:));
    if sum(Pmasktmp(:))>0
        disp(['Displaying result for ',label{ilabel}])
        [maskbw,cluster]=bwlabel(Pmasktmp);
        figure;imagesc(maskbw);axis equal off;colorbar
        for ic = 1:cluster
            disp(['  Numrical report for cluster ' num2str(ic)])
            cluster_sel = maskbw == ic;
            Fvalueall = Fmap(ilabel,cluster_sel);
            
            disp(['    The local maximum within the cluster: F(',num2str(dof(ilabel,1)),...
                ',',num2str(dof(ilabel,2)),') = ',num2str(max(Fvalueall))])
            if isfield(StatMap_c,'beta')
                maxloc  = Fvalueall==max(Fvalueall);
                beta    = betaall(ilabel,cluster_sel);
                betaCI  = squeeze(betaCIall(ilabel,:,cluster_sel));
                disp(['    The beta contrast at the local maximum is ', num2str(beta(maxloc)),...
                    ' [',num2str(betaCI(1,maxloc)),', ',num2str(betaCI(2,maxloc)),']'])
            end
            if plotcondi==1
                for ic2 = 1:size(betamaps,1)
                    beta1=squeeze(betamaps(ic2,:,cluster_sel));
                    maxloc  = Fvalueall==max(Fvalueall);
                    disp(['    Condition ',num2str(ic2),': ',num2str(beta1(1,maxloc)),...
                        ' [',num2str(beta1(2,maxloc)),', ',num2str(beta1(3,maxloc)),']'])
                end
            end
            
            disp(['    The local minimum within the cluster: F(',num2str(dof(ilabel,1)),...
                ',',num2str(dof(ilabel,2)),') = ',num2str(min(Fvalueall))])
            if isfield(StatMap_c,'beta')
                minloc  = Fvalueall==min(Fvalueall);
                beta    = betaall(ilabel,cluster_sel);
                betaCI  = squeeze(betaCIall(ilabel,:,cluster_sel));
                disp(['    The beta contrast at the local minimum is ', num2str(beta(minloc)),...
                    ' [',num2str(betaCI(1,minloc)),', ',num2str(betaCI(2,minloc)),']'])
            end
            if plotcondi==1
                for ic2 = 1:size(betamaps,1)
                    beta1=squeeze(betamaps(ic2,:,cluster_sel));
                    minloc  = Fvalueall==min(Fvalueall);
                    disp(['    Condition ',num2str(ic2),': ',num2str(beta1(1,minloc)),...
                        ' [',num2str(beta1(2,minloc)),', ',num2str(beta1(3,minloc)),']'])
                end
            end
        end
    end
end