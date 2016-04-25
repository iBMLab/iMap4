function normalization(handles,S,choice)

switch choice
    case 1 % duration
        %%
        if (isfield(handles,'FixMap_single_trial_scaled')) == 1 && (isfield(handles,'FixMap_estimated')) == 1
            
            h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
            set(h,'color','white');
            delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
            set(handles.figure1, 'pointer', 'watch')
            drawnow;
            %Single Trial Method
            durind_single_trial = double(handles.duration_single_trial);
            FixMap_duration_single_trial = bsxfun(@rdivide,handles.FixMap_single_trial_scaled,durind_single_trial);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
            FixMap = FixMap_duration_single_trial*100;
            
            
            if isunix ==0
                save(strcat(handles.filename,'\FixMap_normalize_duration_single_trial_scaled'),'FixMap','-v7.3')
                
            else
                save(strcat(handles.filename,'/FixMap_normalize_duration_single_trial_scaled'),'FixMap','-v7.3')
            end
            
            %Estimated Method
            durind_estimated = double(handles.duration_estimated);
            if isfield(handles,'FixMap_estimated_scaled') == 1
                FixMap_duration_estimated = bsxfun(@rdivide,handles.FixMap_estimated_scaled,durind_estimated);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
            else
                FixMap_duration_estimated = bsxfun(@rdivide,handles.FixMap_estimated,durind_estimated);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
            end
            FixMap = FixMap_duration_estimated*100;
            
            
            if isunix ==0
                if isfield(handles,'FixMap_estimated_scaled') == 1
                    save(strcat(handles.filename,'\FixMap_normalize_duration_estimated_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'\FixMap_normalize_duration_estimated'),'FixMap','-v7.3')
                end
            else
                if isfield(handles,'FixMap_estimated_scaled') == 1
                    save(strcat(handles.filename,'/FixMap_normalize_duration_estimated_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'/FixMap_normalize_duration_estimated'),'FixMap','-v7.3')
                end
            end
            set(handles.figure1, 'pointer', 'arrow')
            if ishandle(S.fh)
                delete(S.fh);
            end
            if ishandle(h)
                delete(h)
            end
            
            h = msgbox(strcat('The normalization for the single trial and estimated method is done and saved in',{' '},handles.filename));
            set(h,'color','white');
            
        else if (isfield(handles,'FixMap_single_trial_scaled')) ==1 %make sure the user did the single trial smoothing
                
                %normalization duration
                h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
                set(h,'color','white');
                delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
                set(handles.figure1, 'pointer', 'watch')
                drawnow;
                durind = double(handles.duration_single_trial);
                FixMap_duration_single_trial = bsxfun(@rdivide,handles.FixMap_single_trial_scaled,durind);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
                FixMap = FixMap_duration_single_trial*100;
                %set(handles.figure1, 'pointer', 'arrow')
                
                
                if isunix ==0
                    save(strcat(handles.filename,'\FixMap_normalize_duration_single_trial_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'/FixMap_normalize_duration_single_trial_scaled'),'FixMap','-v7.3')
                    
                end
                set(handles.figure1, 'pointer', 'arrow')
                if ishandle(S.fh)
                    delete(S.fh);
                end
                if ishandle(h)
                    delete(h)
                end
                
                h = msgbox(strcat('The normalization for the single trial method is done and saved in',{' '},handles.filename));
                set(h,'color','white');
                
                
            else  %the user did the estimated smoothing
                %normalization duration
                h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
                set(h,'color','white')
                delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
                set(handles.figure1, 'pointer', 'watch')
                drawnow;
                durind = double(handles.duration_estimated);
                if isfield(handles,'FixMap_estimated_scaled')
                    FixMap_duration_estimated = bsxfun(@rdivide,handles.FixMap_estimated_scaled,durind);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
                else
                    FixMap_duration_estimated = bsxfun(@rdivide,handles.FixMap_estimated,durind);% divided by the number of trials %(smoothpic-mean(smoothpic(:)))./std(smoothpic(:));%
                    
                end
                FixMap = FixMap_duration_estimated*100;
                
                
                if isunix==0
                    if isfield(handles,'FixMap_estimated_scaled') == 1
                        save(strcat(handles.filename,'\FixMap_normalize_duration_estimated_scaled'),'FixMap','-v7.3')
                    else
                        save(strcat(handles.filename,'\FixMap_normalize_duration_estimated'),'FixMap','-v7.3')
                    end
                else
                    if isfield(handles,'FixMap_estimated_scaled') == 1
                        save(strcat(handles.filename,'/FixMap_normalize_duration_estimated_scaled'),'FixMap','-v7.3')
                    else
                        save(strcat(handles.filename,'/FixMap_normalize_duration_estimated'),'FixMap','-v7.3')
                    end
                end
                set(handles.figure1, 'pointer', 'arrow')
                if ishandle(S.fh)
                    delete(S.fh);
                end
                if ishandle(h)
                    delete(h)
                end
                h = msgbox(strcat('The normalization for the estimated method is done and saved in',{' '},handles.filename));
                set(h,'color','white');
            end
        end
        
    case 2 % z_score method
        %%
        if (isfield(handles,'FixMap_single_trial_scaled')) == 1 && (isfield(handles,'FixMap_estimated')) == 1
            
            %Single Trial Method
            FixMap_tmp=handles.FixMap_single_trial_scaled;
            mu = mean(FixMap_tmp(:,:),2);
            sigma = std(FixMap_tmp(:,:),0,2);
            sigma0 = sigma;
            sigma0(sigma0==0) = 1;
            z = bsxfun(@minus,FixMap_tmp, mu);
            FixMap_zscore_single_trial = bsxfun(@rdivide, z,sigma0);
            FixMap = FixMap_zscore_single_trial;
            
            if isunix ==0
                save(strcat(handles.filename,'\FixMap_normalize_zscore_single_trial_scaled'),'FixMap','-v7.3')
                
            else
                save(strcat(handles.filename,'/FixMap_normalize_zscore_single_trial_scaled'),'FixMap','-v7.3')
            end
            
            %Estimated Method
            if isfield(handles,'FixMap_estimated_scaled') == 1
                
                h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
                set(h,'color','white');
                delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
                set(handles.figure1, 'pointer', 'watch')
                drawnow;
                
                FixMap_tmp=handles.FixMap_estimated_scaled;
                mu = mean(FixMap_tmp(:,:),2);
                sigma = std(FixMap_tmp(:,:),0,2);
                sigma0 = sigma;
                sigma0(sigma0==0) = 1;
                z = bsxfun(@minus,FixMap_tmp, mu);
                FixMap_zscore_estimated_scaled = bsxfun(@rdivide, z,sigma0);
                FixMap = FixMap_zscore_estimated_scaled;
            else
                FixMap_tmp=handles.FixMap_estimated;
                mu = mean(FixMap_tmp(:,:),2);
                sigma = std(FixMap_tmp(:,:),0,2);
                sigma0 = sigma;
                sigma0(sigma0==0) = 1;
                z = bsxfun(@minus,FixMap_tmp, mu);
                FixMap_zscore_estimated = bsxfun(@rdivide, z,sigma0);
                FixMap = FixMap_zscore_estimated;
            end
            
            
            if isunix ==0
                if isfield(handles,'FixMap_estimated_scaled') == 1
                    save(strcat(handles.filename,'\FixMap_normalize_zscore_estimated_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'\FixMap_normalize_zscore_estimated'),'FixMap','-v7.3')
                end
            else
                if isfield(handles,'FixMap_estimated_scaled') == 1
                    save(strcat(handles.filename,'/FixMap_normalize_zscore_estimated_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'/FixMap_normalize_zscore_estimated'),'FixMap','-v7.3')
                end
            end
            set(handles.figure1, 'pointer', 'arrow')
            if ishandle(S.fh)
                delete(S.fh);
            end
            if ishandle(h)
                delete(h)
            end
            
            h = msgbox(strcat('The normalization for the single trial and estimated method is done and saved in',{' '},handles.filename));
            set(h,'color','white');
            
        else if (isfield(handles,'FixMap_single_trial_scaled')) ==1 %make sure the user did the single trial smoothing
                h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
                set(h,'color','white');
                delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
                set(handles.figure1, 'pointer', 'watch')
                drawnow;
                %normalization zscore
                FixMap_tmp=handles.FixMap_single_trial_scaled;
                mu = mean(FixMap_tmp(:,:),2);
                sigma = std(FixMap_tmp(:,:),0,2);
                sigma0 = sigma;
                sigma0(sigma0==0) = 1;
                z = bsxfun(@minus,FixMap_tmp, mu);
                FixMap_zscore_single_trial_scaled = bsxfun(@rdivide, z,sigma0);
                FixMap = FixMap_zscore_single_trial_scaled;
                
                
                if isunix ==0
                    save(strcat(handles.filename,'\FixMap_normalize_zscore_single_trial_scaled'),'FixMap','-v7.3')
                else
                    save(strcat(handles.filename,'/FixMap_normalize_zscore_single_trial_scaled'),'FixMap','-v7.3')
                    
                end
                set(handles.figure1, 'pointer', 'arrow')
                if ishandle(S.fh)
                    delete(S.fh);
                end
                if ishandle(h)
                    delete(h)
                end
                
                h = msgbox(strcat('The normalization for the single trial method is done and saved in',{' '},handles.filename));
                set(h,'color','white');
                
                
            else  %the user did the estimated smoothing
                %normalization zscore
                h=msgbox('Please wait... saving the Normalized FixMap might take time ...');
                set(h,'color','white')
                delete(findobj(h,'string','OK')); delete(findobj(h,'style','frame'));
                set(handles.figure1, 'pointer', 'watch')
                drawnow;
                if isfield(handles,'FixMap_estimated_scaled')  == 1
                    FixMap_tmp=handles.FixMap_estimated_scaled;
                    mu = mean(FixMap_tmp(:,:),2);
                    sigma = std(FixMap_tmp(:,:),0,2);
                    sigma0 = sigma;
                    sigma0(sigma0==0) = 1;
                    z = bsxfun(@minus,FixMap_tmp, mu);
                    FixMap_zscore_estimated_scaled = bsxfun(@rdivide, z,sigma0);
                    FixMap = FixMap_zscore_estimated_scaled;
                else
                    FixMap_tmp=handles.FixMap_estimated_scaled;
                    mu = mean(FixMap_tmp(:,:),2);
                    sigma = std(FixMap_tmp(:,:),0,2);
                    sigma0 = sigma;
                    sigma0(sigma0==0) = 1;
                    z = bsxfun(@minus,FixMap_tmp, mu);
                    FixMap_zscore_estimated = bsxfun(@rdivide, z,sigma0);
                    FixMap = FixMap_zscore_estimated;
                    
                end
                
                
                if isunix==0
                    if isfield(handles,'FixMap_estimated_scaled') == 1
                        save(strcat(handles.filename,'\FixMap_normalize_zscore_estimated_scaled'),'FixMap','-v7.3')
                    else
                        save(strcat(handles.filename,'\FixMap_normalize_zscore_estimated'),'FixMap','-v7.3')
                    end
                else
                    if isfield(handles,'FixMap_estimated_scaled') == 1
                        save(strcat(handles.filename,'/FixMap_normalize_zscore_estimated_scaled'),'FixMap','-v7.3')
                    else
                        save(strcat(handles.filename,'/FixMap_normalize_zscore_estimated'),'FixMap','-v7.3')
                    end
                end
                set(handles.figure1, 'pointer', 'arrow')
                if ishandle(S.fh)
                    delete(S.fh);
                end
                if ishandle(h)
                    delete(h)
                end
                h = msgbox(strcat('The normalization for the estimated method is done and saved in',{' '},handles.filename));
                set(h,'color','white');
            end
        end
        
        
        
end