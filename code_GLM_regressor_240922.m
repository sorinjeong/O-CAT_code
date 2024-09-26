clear; clc; close;
path='G:\JSR\240922_new_fmri\bhv\data_bhv_log_table\total';
addpath(path)
 
root_path = 'G:\JSR'; %69번컴
% addpath(genpath(root_path))
cd(fullfile(root_path, 'sbj39to41'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%% make events log table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('regressors_GLM_0922_backup.mat')
% %% edit event table --> add Navigation period
% load('G:\JSR\sbj39to41\sbj_id_list_38.mat');
% subjects=[sbj_id_list_38 39:41];
% files=dir(fullfile(path, '*events.tsv'));
% for i=1:numel(files)
% T=readtable(files(i).name,"FileType","text",'Delimiter', '\t');
% T.TrialStart(T.Lap~=0 & T.Lap~=9)=T.TrialStart(T.Lap~=0 & T.Lap~=9)-4;
% T=renamevars(T,["TrialStart", "TrialEnd"],["NavStart", "ITIEnd"]);
% T.NavEnd=NaN(height(T),1);
% T.NavEnd(T.Lap_Trial==4)=T.ITIEnd(T.Lap_Trial==4)+4;
% writetable(T,'all_subjects_task-ocat_run-01_events.xlsx','Sheet',sprintf('sub-%.2d',i))
% end
% 
% 
% % align TR with events
% x="all_subjects_task-ocat_run-01_events.xlsx";
% 
% [~,sheets] = xlsfinfo(x);
% 
% for s=1:numel(sheets)
%     table_events= readtable(x,'Sheet',sheets{s});
%     table_TR= readtable('event_TR.xlsx','Sheet',sheets{s});
% 
%     table_events.MR_start=zeros(height(table_events),1);
%     table_events.MR_end=zeros(height(table_events),1);
% 
%     for i= 1:height(table_events)
%         [~,idx1]=min(abs(table_TR.Var2 - table_events.NavStart(i)));
%         if table_events.Trial(i) <0
%             [~,idx2]=min(abs(table_TR.Var2 - table_events.ObjOff(i)));
%         else
%             [~,idx2]=min(abs(table_TR.Var2 - table_events.NavEnd(i)));
%         end
%         table_events.MR_start(i)=table_TR.Var2(idx1);
%         table_events.MR_end(i)=table_TR.Var2(idx2);
% 
%     end
% 
%     writetable(table_events, fullfile(path, 'event_table_MR.xlsx'),'Sheet',sheets{s});
% end
% 
% % make TR scan
% cd('G:\JSR\240922_new_fmri\bhv\data_bhv_log_table\total')
% %%
% sheet_names=sheetnames("event_TR.xlsx");
% numsheets=numel(sheet_names);
% for i= 1: numsheets
% c_sbj=sprintf('sub-%02d',i);
% curr_T= readtable("event_TR.xlsx","Sheet",c_sbj);
% %%
% first_scan=find(curr_T.Var2==0,1,'last');
% curr_T.scan_num=[zeros(first_scan-1,1);(1:height(curr_T)-first_scan+1)'];
% writetable(curr_T,"event_TR_scan.xlsx","Sheet",c_sbj)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %% add ODT object off onset
% for n_sub = 1:41
%     c_sub=sprintf('sub-%.2d',n_sub);
% 
% A = readtable('event_pre_ODT.xlsx', 'Sheet', c_sub);
% B = readtable('event_table_MR.xlsx', 'Sheet', c_sub);
% C = readtable('event_post_ODT.xlsx', 'Sheet', c_sub);
% 
% 
% preoff = A.Var2(strcmp(A.Var1,'PreObjOff'));
% B.ObjOff(1:15) = preoff;
% postoff = C.Var2(strcmp(C.Var1,'PreObjOff'));
% B.ObjOff(48:62) = postoff;
% 
% writetable(B, 'event_table_MR.xlsx', 'Sheet', c_sub);
% end
%%%%%%% make regressor cell array %%%%%%%%%%%%%%%%%%
% sbj_id_list_41 = subjects;
%  load('scans_38.mat')
%  scans=[scans 375 375 375];
% save("scans_41","scans");
% scans_37=scans(setdiff(1:41,[2,3,7,21]));
% save("sbj_id_list_41","sbj_id_list_41", "scans_37");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%make regressor struct%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT!!
TR = 2;

load('sbj_id_list_41');
scans=scans_37;

%% OD task: regressor names cell
ODT_regress_name=cell(2,5);ODT_regress_duration=cell(2,5);
% ODT=struct('pre',struct,'post',struct);
%name
ODT_regress_name(1,:) = {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','pre_target'};
ODT_regress_name(2,:) = {'post_forest_1','post_forest_2','post_city_1','post_city_2','post_target'};

%duration
% pre-ODT_all = 120s; post-ODT_all = 120s
ODT_regress_duration(1,:) = repmat({repelem(4,1,3)},1,5);
ODT_regress_duration(2,:) = repmat({repelem(4,1,3)},1,5);

%%
reg_0922={};ODT_regress_onset=cell(2,5);ODT=struct;
for i = 1: numel(sbj_id_list_41)
    n_sbj = sbj_id_list_41(i);
    c_sbj = sprintf('sub-%.2d',n_sbj); disp(c_sbj)
    
    curr_T = readtable([path '\event_table_MR.xlsx'],'Sheet',c_sbj);
    curr_TR = readtable([path '\event_TR_scan.xlsx'],'Sheet',c_sbj);
    
    % onset after trim
    main_onset= curr_T.NavStart(curr_T.Trial==1)-6;
    main_offset = curr_T.NavEnd(curr_T.Trial==32);
    %     pre_onset = curr_T.ObjOn(curr_T.Trial==-1);
    pre_offset = curr_T.ObjOff(curr_T.Trial==-15);
    post_onset = curr_T.ObjOn(curr_T.Trial==-16);
    %     post_offset = curr_T.ObjOff(curr_T.Trial==-30);
    
    % scan
    main_scan_num=[floor(main_onset/TR), round(main_offset/TR)+3];
    %     pre_scan_num=[floor(pre_onset/TR), round(pre_offset/TR)+3];
    %     post_scan_num=[floor(post_onset/TR), round(post_offset/TR)+3];
    pre_scan_num=[1, round(pre_offset/TR)+5];
    post_scan_num=[round(post_onset/TR)-5, scans(i)];
    
    % scan onset time
    % main_scan_onset=curr_TR.Var2(floor(main_onset/TR)+2);
    % pre_scan_onset=curr_TR.Var2(floor(pre_onset/TR)+2);
    % post_scan_onset=curr_TR.Var2(floor(post_onset/TR)+2);
    
    starting_scan = find(curr_TR.Var2==0,1,"last");
    % onset time change
    main_scan_onset = curr_TR.Var2(floor(main_onset/TR)+starting_scan-1); % fourth TR
    %     pre_scan_onset = curr_TR.Var2(floor(pre_onset/TR)+starting_scan-1);
    pre_scan_onset = curr_TR.Var2(curr_TR.scan_num==1);
    
    post_scan_onset = curr_TR.Var2((round(post_onset/TR)-5)+starting_scan-1);
    if ~((0 < abs(main_scan_onset - (floor(main_onset/TR)-1)*TR)) && (0.5 > abs(main_scan_onset - (floor(main_onset/TR)-1)*TR)))
        error('scan num error');
    end
    
    %% onset
    % associated object
    obj_f = unique(curr_T.Obj_ID(curr_T.Context_Num ==1 & curr_T.Association == 1));
    obj_c = unique(curr_T.Obj_ID(curr_T.Context_Num ==2 & curr_T.Association == 1));
    
    f1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_f) & curr_T.Trial < 0);
    f2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_f) & curr_T.Trial < 0);
    c1 = curr_T.ObjOn(curr_T.Obj_ID == min(obj_c) & curr_T.Trial < 0);
    c2 = curr_T.ObjOn(curr_T.Obj_ID == max(obj_c) & curr_T.Trial < 0);
    
    % ODT
    for odt_row=1:2
        if odt_row == 1; obj_idx=1:3; ons=pre_scan_onset;
        elseif odt_row == 2; obj_idx = 4:6; ons=post_scan_onset;
        end
        ODT_regress_onset{odt_row,1} =f1(obj_idx)-ons;
        ODT_regress_onset{odt_row,2} =f2(obj_idx)-ons;
        ODT_regress_onset{odt_row,3} =c1(obj_idx)-ons;
        ODT_regress_onset{odt_row,4} =c2(obj_idx)-ons;
    end
    target = curr_T.ObjOn(curr_T.Obj_ID==12);
    ODT_regress_onset{1,5} = target(1:3) - pre_scan_onset;
    ODT_regress_onset{2,5} = target(4:6)  - post_scan_onset;
    
    
    
    %% Main
    main_T=curr_T(~isnan(curr_T.Correct_Num),:);
    onset=struct;
    
    % time out
    idx_to = (main_T.Correct_Num == 2);
    if any(idx_to); trial_detail.to=find(idx_to); end
    
    %% object cueing, choice ROC
    fix_roc={};
    for cols={'ObjOn' 'ChoiceOn'}
        % Hit
        idx_hit = (main_T.Correct_Num == 1) & (main_T.Association == 1);
        if any(idx_hit); onset.hit.(cols{:}) = table2array(main_T(idx_hit, cols{:}));
            fix_roc{1}='hit';
            % trial_detail.hit=setdiff(find(idx_hit), find(idx_to)); end
            trial_detail.hit=find(idx_hit); end
        
        
        % correct rejection
        idx_corr = (main_T.Correct_Num == 1) & (main_T.Association == 0);
        if any(idx_corr); onset.corr_rej.(cols{:}) = table2array(main_T(idx_corr, cols{:}));
            fix_roc{4}='corr_rej';
            % trial_detail.corr_rej=setdiff(find(idx_corr), find(idx_to)); end
            trial_detail.corr_rej=find(idx_corr); end
        
        
        % miss
        idx_miss = (main_T.Correct_Num ~= 1) & (main_T.Association == 1);
        if any(idx_miss); onset.miss.(cols{:}) = table2array(main_T(idx_miss, cols{:}));
            fix_roc{2}='miss';
            % trial_detail.miss=setdiff(find(idx_miss), find(idx_to)); end
            trial_detail.miss=find(idx_miss); end
        
        
        % false positive
        idx_false = (main_T.Correct_Num ~= 1) & (main_T.Association == 0);
        if any(idx_false); onset.false.(cols{:}) = table2array(main_T(idx_false, cols{:}));
            fix_roc{3}='false';
            % trial_detail.false=setdiff(find(idx_false), find(idx_to)); end
            trial_detail.false=find(idx_false); end
        
    end
    
    %% choice - button press
    idx_L = (main_T.Choice_Num == 1);
    idx_R = (main_T.Choice_Num == 2);
    
    if any(idx_L); onset.button.A = table2array(main_T(idx_L, 'ChoiceOn'));end
    if any(idx_R); onset.button.B = table2array(main_T(idx_R, 'ChoiceOn'));end
    
    
    %% Moving phase
    
    idx_mov_1_f = (main_T.Lap_Trial == 1 & main_T.Context_Num == 1);
    idx_mov_3_f = (main_T.Lap_Trial == 3 & main_T.Context_Num == 1);
    
    if any(idx_mov_1_f); onset.mov.f.loc_1 = main_T.NavStart(idx_mov_1_f)+2;end
    if any(idx_mov_1_f); onset.mov.f.loc_2 = main_T.ITIEnd(idx_mov_1_f)+2;end
    if any(idx_mov_3_f); onset.mov.f.loc_3 = main_T.ITIEnd(idx_mov_3_f)+2;end
    if any(idx_mov_3_f); onset.mov.f.loc_4 = main_T.NavStart(idx_mov_3_f)+2;end
    
    idx_mov_1_c = (main_T.Lap_Trial == 1 & main_T.Context_Num == 2);
    idx_mov_3_c = (main_T.Lap_Trial == 3 & main_T.Context_Num == 2);
    
    if any(idx_mov_1_c); onset.mov.c.loc_1 = main_T.NavStart(idx_mov_1_c)+2;end
    if any(idx_mov_1_c); onset.mov.c.loc_2 = main_T.ITIEnd(idx_mov_1_c)+2;end
    if any(idx_mov_3_c); onset.mov.c.loc_3 = main_T.ITIEnd(idx_mov_3_c)+2;end
    if any(idx_mov_3_c); onset.mov.c.loc_4 = main_T.NavStart(idx_mov_3_c)+2;end
    
    %% main task: regressor names cell
    main_regress_name={};main_regress_duration={};
    %% regressor name
    main_regress_name = {'obj_show_hit','obj_show_miss','obj_show_false','obj_show_corr_rej','choice_A','choice_B','mov'}; %new best
    
    %% regressor onset
    fix_phase={'ObjOn' 'ChoiceOn','mov'};
    
    fix_roc(cellfun(@(x) isnumeric(x) && isempty(x), fix_roc))=[];
    oc_onset = cellfun(@(x) onset.(x).(fix_phase{1}), fix_roc, 'UniformOutput', false);
    oc_onset = cellfun(@(x) x(~isnan(x))-main_scan_onset, oc_onset, 'UniformOutput', false);
    
    ch_onset = cellfun(@(x) onset.(x).(fix_phase{2}), fix_roc, 'UniformOutput', false);
    ch_onset = cellfun(@(x) x(~isnan(x))-main_scan_onset, ch_onset, 'UniformOutput', false);
    
    button_onset = {onset.button.A onset.button.B};
    button_onset = cellfun(@(x) x(~isnan(x))-main_scan_onset, button_onset, 'UniformOutput', false);
    
    mov_onset = horzcat(arrayfun(@(x)  onset.(fix_phase{3}).f.(['loc_' num2str(x)]), 1:4, 'UniformOutput', false), arrayfun(@(x)  onset.(fix_phase{3}).c.(['loc_' num2str(x)]), 1:4, 'UniformOutput', false));
    mov_onset = cellfun(@(x) x(~isnan(x))-main_scan_onset, mov_onset, 'UniformOutput', false);
    
    
    main_regress_onset = [oc_onset button_onset {sort(vertcat(mov_onset{:})')}];


        main.name=main_regress_name;
        main.onset=main_regress_onset;
        main.duration=zeros(size(main_regress_onset));
        main.scan_num=main_scan_num;

    
    
    %% set onset to scan start time
    
    ODT.pre= struct('scan_num',pre_scan_num,'regress_name',{ODT_regress_name(1,:)'},...
        'regress_onset',{ODT_regress_onset(1,:)'},'regress_duration',{ODT_regress_duration(1,:)'});
    ODT.post= struct('scan_num',post_scan_num,'regress_name',{ODT_regress_name(2,:)'},...
        'regress_onset',{ODT_regress_onset(2,:)'},'regress_duration',{ODT_regress_duration(2,:)'});
    
    
    reg_0922{1,i} = struct('n_sbj',n_sbj,'main',main,'ODT',ODT,'trial_detail',trial_detail);
    
end

save(fullfile(root_path,'240922_new_fmri','regressors_GLM_0922.mat'),"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name');
save('regressors_GLM_0922_backup.mat',"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name');
