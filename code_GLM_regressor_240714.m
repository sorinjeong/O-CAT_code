clear; clc; close;
path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_BHV\data\0716\data_bhv_log_table_38\total';
cd(path); load('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code\sbj_id_list.mat');

%%%%%%%% make events log table %%%%%%%%%%

%% edit event table --> add Navigation period
files=dir(fullfile(path, '*events.tsv'));
for i=1:numel(files)
T=readtable(files(i).name,"FileType","text",'Delimiter', '\t');
T.TrialStart(T.Lap~=0 & T.Lap~=9)=T.TrialStart(T.Lap~=0 & T.Lap~=9)-4;
T=renamevars(T,["TrialStart", "TrialEnd"],["NavStart", "ITIEnd"]);
T.NavEnd=NaN(height(T),1);
T.NavEnd(T.Lap_Trial==4)=T.ITIEnd(T.Lap_Trial==4)+4;
writetable(T,'all_subjects_task-ocat_run-01_events.xlsx','Sheet',sprintf('sub-%.2d',i))
end


% align TR with events
x="all_subjects_task-ocat_run-01_events.xlsx";

[~,sheets] = xlsfinfo(x);

for s=1:numel(sheets)
    table_events= readtable(x,'Sheet',sheets{s});
    table_TR= readtable('event_TR.xlsx','Sheet',sheets{s});

    table_events.MR_start=zeros(height(table_events),1);
    table_events.MR_end=zeros(height(table_events),1);

    for i= 1:height(table_events)
        [~,idx1]=min(abs(table_TR.Var2 - table_events.NavStart(i)));
        if table_events.Trial(i) <0
            [~,idx2]=min(abs(table_TR.Var2 - table_events.ObjOff(i)));
        else
            [~,idx2]=min(abs(table_TR.Var2 - table_events.NavEnd(i)));
        end
        table_events.MR_start(i)=table_TR.Var2(idx1);
        table_events.MR_end(i)=table_TR.Var2(idx2);

    end

    writetable(table_events, 'event_table_MR.xlsx','Sheet',sheets{s});
end
%
%% add ODT object off onset
for n_sub = 1:38
    c_sub=sprintf('sub-%.2d',n_sub);

A = readtable('event_pre_PV.xlsx', 'Sheet', c_sub);
B = readtable('event_table_MR.xlsx', 'Sheet', c_sub);
C = readtable('event_post_PV.xlsx', 'Sheet', c_sub);


preoff = A.Var2(strcmp(A.Var1,'PreObjOff'));
B.ObjOff(1:15) = preoff;
postoff = C.Var2(strcmp(C.Var1,'PreObjOff'));
B.ObjOff(48:62) = postoff;

writetable(B, 'event_table_MR.xlsx', 'Sheet', c_sub);
end

%%%%%%% make regressor cell array %%%%%%%%%%%%%%%%%%
sbj_id_list_38 = [sbj_id_list 32:38];
save("sbj_id_list_38","sbj_id_list_38")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%make regressor struct%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT!!
TR = 2;

load('sbj_id_list_38');
%% main task: regressor names cell
main_regress_name={};main_regress_duration={};
% name
main_regress_name(1,:) = {'obj_show_hit','obj_show_miss','obj_show_false','obj_show_corr_rej','obj_ITI'};
main_regress_name(2,:) = {'choice_hit','choice_miss','choice_false','choice_corr_rej','moving'};
% duration
main_regress_duration(1,:) = {4,4,4,4,2};
main_regress_duration(2,:) = {2,2,2,2,4};

%% OD task: regressor names cell
ODT_regress_name=cell(2,5);ODT_regress_duration=cell(2,5);
% ODT=struct('pre',struct,'post',struct);
%name
ODT_regress_name(1,:) = {'pre_forest_1','pre_forest_2','pre_city_1','pre_city_2','pre_target'};
ODT_regress_name(2,:) = {'post_forest_1','post_forest_2','post_city_1','post_city_2','post_target'};

%duration
% pre-ODT_all = 120s; post-ODT_all = 120s
ODT_regress_duration(1,:) = {4,4,4,4,4};
ODT_regress_duration(2,:) = {4,4,4,4,4};

%%
reg_0717={};ODT_regress_onset=cell(2,5);ODT=struct;
for i = 1: numel(sbj_id_list_38)
    n_sbj = sbj_id_list_38(i);
    c_sbj = sprintf('sub-%.2d',n_sbj);

    curr_T = readtable('event_table_MR.xlsx','Sheet',c_sbj);

    curr_TR = readtable('event_TR.xlsx','Sheet',c_sbj);

    % onset after trim
    main_onset= curr_T.NavStart(curr_T.Trial==1);
    main_offset = curr_T.NavEnd(curr_T.Trial==32);
    pre_onset = curr_T.ObjOn(curr_T.Trial==-1);
    pre_offset = curr_T.ObjOff(curr_T.Trial==-15);
    post_onset = curr_T.ObjOn(curr_T.Trial==-16);
    post_offset = curr_T.ObjOff(curr_T.Trial==-30);

    % scan
    
    main_scan_num=[floor(main_onset/TR), round(main_offset/TR)+2];
    pre_scan_num=[floor(pre_onset/TR), round(pre_offset/TR)+2];
    post_scan_num=[floor(post_onset/TR), round(post_offset/TR)+2];

% scan onset time
% main_scan_onset=curr_TR.Var2(floor(main_onset/TR)+2);
% pre_scan_onset=curr_TR.Var2(floor(pre_onset/TR)+2);
% post_scan_onset=curr_TR.Var2(floor(post_onset/TR)+2);

    starting_scan = find(curr_TR.Var2==0,1,"last");
    % onset time change
    main_scan_onset = curr_TR.Var2(floor(main_onset/TR)+starting_scan-1); % fourth TR
    pre_scan_onset = curr_TR.Var2(floor(pre_onset/TR)+starting_scan-1);
    post_scan_onset = curr_TR.Var2(floor(post_onset/TR)+starting_scan-1);
    if ~((0 < abs(main_scan_onset - (floor(main_onset/TR)-1)*TR)) && (0.5 > abs(main_scan_onset - (floor(main_onset/TR)-1)*TR)))
        break;
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


    % onset
    onsets=[];main_regress_onset=cell(2,5);fix=[];
    trial_hit=[];trial_corr_rej=[];trial_miss=[];trial_false=[];trial_to=[];


    for row = 1:height(curr_T)
        %% main task
        for reg_row = 1:2
            if reg_row == 1; onsets = curr_T.ObjOn;
            elseif reg_row == 2; onsets = curr_T.ChoiceOn;
            end
if ~isnan(curr_T.Correct_Num(row))
            if curr_T.Correct_Num(row) ==1
                % Hit
                if curr_T.Association(row) ==1
                    main_regress_onset{reg_row,1}(end+1) = onsets(row);
                    if reg_row == 1;trial_hit(end+1) = curr_T.Trial(row);end
                    % Correct_rejection
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,4}(end+1) = onsets(row);
                    if reg_row == 1;trial_corr_rej(end+1) = curr_T.Trial(row);end

                end
            else
                % miss
                if curr_T.Association(row) ==1
                    main_regress_onset{reg_row,2}(end+1) = onsets(row);
                    if reg_row == 1;trial_miss(end+1) = curr_T.Trial(row);end

                    % false alarm
                elseif curr_T.Association(row) ==0
                    main_regress_onset{reg_row,3}(end+1) = onsets(row);
                    if reg_row == 1;trial_false(end+1) = curr_T.Trial(row);end

                end

            end
end
        end
        if curr_T.Trial(row) > 0

            % Moving
            if curr_T.Lap_Trial(row) ==1
                main_regress_onset{2,5}(end+1) = curr_T.NavStart(row);
            end

            % object ITI
            main_regress_onset{1,5}(end+1) = curr_T.ObjOff(row);
   
        end


    end

    %% set onset to scan start time

    main_regress_onset = cellfun(@(x) x(~isnan(x))-main_scan_onset, main_regress_onset, 'UniformOutput', false);
names = reshape(main_regress_name',[],1);
onsets = reshape(main_regress_onset',[],1);
durations = reshape(main_regress_duration',[],1);
order = [1,2,3,4,6,7,8,9,5,10];
    %%

    main=struct('scan_num',main_scan_num,'regress_name',{names(order)},...
        'regress_onset',{onsets(order)},'regress_duration',{durations(order)});
     % ODT = struct('pre_scan_num',pre_scan_num,'post_scan_num',post_scan_num,'regress_name',{ODT_regress_name},...
     %    'regress_onset',{ODT_regress_onset},'regress_duration',{ODT_regress_duration});


    ODT.pre= struct('scan_num',pre_scan_num,'regress_name',{ODT_regress_name(1,:)'},...
        'regress_onset',{ODT_regress_onset(1,:)'},'regress_duration',{ODT_regress_duration(1,:)'});
    ODT.post= struct('scan_num',post_scan_num,'regress_name',{ODT_regress_name(2,:)'},...
        'regress_onset',{ODT_regress_onset(2,:)'},'regress_duration',{ODT_regress_duration(2,:)'});
    

    reg_0717{1,i} = struct('n_sbj',n_sbj,'main',main,'ODT',ODT,'trial_detail',...
        struct('trial_hit',trial_hit,'trial_corr_rej',trial_corr_rej,'trial_miss',trial_miss,'trial_false',trial_false,'trial_to',trial_to));

end
root_path = 'Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
save(fullfile(root_path,'regressors_GLM_0717'),"reg_0717","sbj_id_list_38");
