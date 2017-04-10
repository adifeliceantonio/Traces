%A. DiFeliceantonio March 10, 2017
%This is a function to analyze extracted traces from calcium imaging ROIs.
%To run this function you need an excel sheet saved as comma deliminated
%.csv that cointains: labels in the first row, time stamps in the first
%column and ROI traces in all subsequent columns. 
%All files should be in the same folder as this script and be named
%RAT_acq.csv and RAT_NCD.csv
%%
%Get information from the experimenter
clear
%Get the rat number
rat_id = inputdlg('Rat Number?', 'Rat Number', 1);
rat_id=rat_id{1,1};
%Get NCD information
ncd_prompts={'Last time frame of NL session', 'Last time frame of post NL session', 'Last time from of 6hrs post NL session'};
ncd_dlg_title='NCD session ifo';
numlines=1;
nl_info=inputdlg(ncd_prompts, ncd_dlg_title, numlines);
nl = str2double(nl_info{1,1});
nl_post=str2double(nl_info{2,1});
nl_6hrs=str2double(nl_info{3,1});

%Get ACQ information
acq_prompts={'Last time frame of Acquisition session', 'Last time frame of post Acquisition session', 'Last time from of 6hrs post Acquisition session'};
acq_dlg_title='Acquisition session ifo';
acq_info=inputdlg(acq_prompts, acq_dlg_title, numlines);
acq = str2double(acq_info{1,1});
acq_post=str2double(acq_info{2,1});
acq_6hrs=str2double(acq_info{3,1});
% Get the data files
acqfilename=sprintf('%s_ACQ.csv', rat_id);
nlfilename=sprintf('%s_NL.csv', rat_id);
acq_session=csvread(acqfilename, 1, 0);
nl_session=csvread(nlfilename, 1,0);

%% Create Data Arrays and Block Regressors
%determine number of arrays and time points
[acq_totaltime, acq_totalrois]=size(acq_session);
[nl_totaltime, nl_totalrois]=size(nl_session);

%create data table for acq session
acq_variablenames=string(1:acq_totalrois-1);
acq_variablenames=strcat('roi_', acq_variablenames);
acq_variablenames=cellstr(acq_variablenames);
acq_data=array2table(acq_session(:,2:acq_totalrois), 'VariableNames', acq_variablenames);
%Create table to write results into for acquisistion session
acq_columnames={'ACQ' 'ACQ_post'};
acq_results=array2table(zeros(acq_totalrois-1, 2), 'RowNames', acq_variablenames, 'VariableNames', acq_columnames);

%create data table for nl sesion
nl_variablenames=string(1:nl_totalrois-1);
nl_variablenames=strcat('roi_', nl_variablenames);
nl_variablenames=cellstr(nl_variablenames);
nl_data=array2table(nl_session(:,2:nl_totalrois), 'VariableNames', nl_variablenames);
%Create table to write results into for acquisistion session
nl_columnames={'NL' 'NL_post'};
nl_results=array2table(zeros(nl_totalrois-1, 2), 'RowNames', nl_variablenames, 'VariableNames', nl_columnames);


%create block regressors
%Acquisistion session
acq_data.blockacq=zeros(acq_totaltime, 1);
findacq=find(acq_session(:,1)<acq);
acq_data.blockacq(findacq,:)=1;
acq_data.blockacq=nominal(acq_data.blockacq);
%10 minutes post acquisistion
acq_data.blockacq_post=zeros(acq_totaltime, 1);
findacq_post=find(acq_session(:,1)>acq & acq_session(:,1)<acq_post);
acq_data.blockacq_post(findacq_post,:)=1;
acq_data.blockacq_post=nominal(acq_data.blockacq_post);
%6hrs post acquisition
acq_data.blockacq_6hrs=zeros(acq_totaltime, 1);
findacq_6hrs=find(acq_session(:,1)>acq_post & acq_session(:,1)<acq_6hrs);
acq_data.blockacq_6hrs(findacq_6hrs,:)=1;
acq_data.blockacq_6hrs=nominal(acq_data.blockacq_6hrs);

%NL
nl_data.blocknl=zeros(nl_totaltime, 1);
findnl=find(nl_session(:,1)<=nl);
nl_data.blocknl(findnl,:)=1;
nl_data.blocknl=nominal(nl_data.blocknl);
%10 minutes post NL 
nl_data.blocknl_post=zeros(nl_totaltime, 1);
findnl_post=find(nl_session(:,1)>nl & nl_session(:,1)<nl_post);
nl_data.blocknl_post(findnl_post,:)=1;
nl_data.blocknl_post=nominal(nl_data.blocknl_post);
%6 hours post NL 
nl_data.blocknl_6hrs=zeros(nl_totaltime, 1);
findmag_6hrs=find(nl_session(:,1)>nl_post & nl_session(:,1)<nl_6hrs);
nl_data.blocknl_6hrs(findmag_6hrs,:)=1;
nl_data.blocknl_6hrs=nominal(nl_data.blocknl_6hrs);


%% Let's DO STATISTICS
%Massively univariate approach to identify ROIS that differ from baseline
%and from the different sessions
H=[0 1 -1];
nl_contrastnames={'acq', 'post_acq'};
nl_finders={findnl, findnl_post};
%run analysis for each ROI in NL session
for i=1:nl_totalrois-1
    roi=sprintf('roi_%d', i);
    nl_model=sprintf('%s~blocknl+blocknl_post', roi);
    nl_lm=fitlm(nl_data, nl_model);  
    nl_modelname=sprintf('nl_%s', roi);
    nl_lms.(nl_modelname)=nl_lm;
    %nl_anovas.(nl_modelname)=anova(nl_lms.(nl_modelname));
    for contrast=1:2
        roi_means.(nl_contrastnames{1,contrast})(i,1)=mean(nl_data.(roi)(nl_finders{contrast},1));
        nl_pvalue.(nl_modelname).(nl_contrastnames{1,contrast})=coefTest(nl_lms.(nl_modelname), H);
    end
    if nl_pvalue.(nl_modelname).(nl_contrastnames{1,contrast})<0.0001
        if roi_means.(nl_contrastnames{1,1})(i,1)> roi_means.(nl_contrastnames{1,2})(i,1)
        nl_results.(nl_columnames{1})(i,1)=1;
        else
        nl_results.(nl_columnames{2})(i,1)=1;
        end      
    end
end
%Massively univariate approach to identify ROIS that differ from baseline
%and from the different sessions
acq_contrastnames={'acq', 'acq_post'};
acq_finders={findacq, findacq_post};
%run analysis for each ROI in acq session
for i=1:acq_totalrois-1
    roi=sprintf('roi_%d', i);
    acq_model=sprintf('%s~blockacq+blockacq_post', roi);
    acq_lm=fitlm(acq_data, acq_model);  
    acq_modelname=sprintf('acq_%s', roi);
    acq_lms.(acq_modelname)=acq_lm;
    %acq_anovas.(acq_modelname)=anova(acq_lms.(acq_modelname));
    for contrast=1:2
        roi_means.(acq_contrastnames{1,contrast})(i,1)=mean(acq_data.(roi)(acq_finders{contrast},1));
        acq_pvalue.(acq_modelname).(acq_contrastnames{1,contrast})=coefTest(acq_lms.(acq_modelname), H);
    end
    if acq_pvalue.(acq_modelname).(acq_contrastnames{1,contrast})<0.0001
        if roi_means.(acq_contrastnames{1,1})(i,1)> roi_means.(acq_contrastnames{1,2})(i,1)
        acq_results.(acq_columnames{1})(i,1)=1;
        else
        acq_results.(acq_columnames{2})(i,1)=1;
        end      
    end
end



% Make Pie charts
labels={'inactive', 'active'};
subplot(2,2,1);
resacq=nominal(acq_results.ACQ);
pie(resacq)
title('ACQ')
subplot(2,2,3);
resacq_post=nominal(acq_results.ACQ_post);
pie(resacq_post)
title('Post ACQ')
subplot(2,2,2)
resnl=nominal(nl_results.NL);
pie(resnl)
title('NL')
subplot(2,2,4)
resnl_post=nominal(nl_results.NL_post);
pie(resnl_post)
title('Post NL')

save_filenameacq=sprintf('%s_ACQ_results.mat', rat_id);
save(save_filenameacq, 'acq_results', 'acq*', '-v7.3');
save_filenamenl=sprintf('%s_NL_results.mat', rat_id);
save(save_filenamenl, 'nl*', '-v7.3');
