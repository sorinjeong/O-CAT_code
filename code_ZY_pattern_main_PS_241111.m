%% single trial naming
clc;clear

% zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
% addpath(zcode_path)
% root_path = 'C:\Users\User\Desktop\JSR';%90번컴
cd('G:\JSR\241104_new_fmri\pattern_similarity');addpath(genpath('G:\JSR\241104_new_fmri'));
root_path = 'G:\JSR\241104_new_fmri'; %69번컴
addpath(fullfile(root_path, 'glm_new_1111')); %addpath('G:\JSR\sbj39to41')
% cd(fullfile(root_path, 'spm_prep_glm_0922','main'))
% zroot_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\data\spm_prep_glm_0725';
% addpath(zroot_path)

load(fullfile(root_path,'bhv/num_sbj_events.mat'))
load('reg_1111.mat')
load(fullfile(root_path,'bhv/learn_curv_info.mat'));


tbl=num_sbj_events;
% setting
main_or_ODT = {'main'}; % {'main','pre-ODT','post-ODT'} for for loops


hpc_select=[1:5,7:17,19:24];
ctx_select=1:30;


%
%%
%
% for m_o = main_or_ODT
%     disp(m_o)
%
%         singletrial_path = fullfile(root_path, 'spm_prep_glm_0922', m_o{:}, 'single_trial');
%
%     %     for sbj_i=1:numel(sbj_id_list)
%     %                 sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
%     %     files=dir(fullfile(singletrial_path,sbj_n, '*_0*'));
%     %     for i=1:height(files)
%     %     delete(fullfile(files(1).folder,files(i).name))
%     %     end
%     %     end
%     %     end
%

singletrial_path=fullfile(root_path, 'glm_new_1111','main','single_trial');

for sbj_i=1:numel(sbj_id_list)
    sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
    
    
    pp='obj_';
    %     sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001',strcat(pp,'*_*')));
    %11/15 수정 윗줄 코드 실수;;;
        sbj_OC=dir(fullfile(singletrial_path,sbj_n,'*.nii'));

    roc=cellfun(@(x) extractBetween(x,'obj_','_0'),{sbj_OC.name});
    d=reg_1111{1, sbj_i}.trial_detail;
    un=unique(roc);
    
    sidx=[];
    for r=1:numel(un)
        idx=d.(un{r})';
        sidx=[sidx idx];
    end
    [b,iii]=sort(sidx);
    if max(iii)==32
        tbl.roc(tbl.session==sbj_id_list(sbj_i))=roc(iii);
        
        %           for i=1:height(sbj_OC)
        %                     beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
        %                     beta_rename=sprintf('%.2d_%s.nii',sidx(i), sbj_OC(i).name);
        %                     copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));
        %           end
        
    else
        error('The total number of trials is not 32!');
    end
end


%% Trial detail
% load("sbj_events.mat")
m_o = main_or_ODT;

disp('trial detail')

roc_cell=cell(1, ceil(length(tbl.roc)/32));
for i = 1:length(roc_cell)
    startIdx = (i-1)*32 + 1;
    endIdx = min(i*32, length(tbl.roc));
    roc_cell{i} = tbl.roc(startIdx:endIdx);
end

% trial별 object 및 context 정보
obj_num=cell(1, numel(sbj_id_list));
context=cell(1, numel(sbj_id_list));
for sbj_i=1:numel(sbj_id_list)
    
    %             startIdx = (sbj_id_list(sbj_i)-1)*32 + 1;
    %             endIdx = min(sbj_id_list(sbj_i)*32, height(tbl));
    obj_num{sbj_i} = tbl.Object(tbl.session==sbj_id_list(sbj_i))';
    t_ctxt=tbl.Context(tbl.session==sbj_id_list(sbj_i))';
    t_asso=tbl.Association(tbl.session==sbj_id_list(sbj_i))';
    for cont=1:32
        if t_asso(cont) ==1
            ctxt{cont} = t_ctxt(cont);
        elseif t_asso(cont) ==0 && t_ctxt(cont) ==1
            ctxt{cont} = 2;
        elseif t_asso(cont) ==0 && t_ctxt(cont) ==2
            ctxt{cont} = 1;
        end
    end
    context{sbj_i}=ctxt;
    
    reg_1111{1, sbj_i}.trial_detail.context=cell2mat(context{sbj_i});
    reg_1111{1, sbj_i}.trial_detail.object=obj_num{sbj_i};
end
save('regressor_1111.mat',"reg_1111","sbj_id_list",'roi_ctx_name','roi_hpc_name')




%% pattern structure
disp('pattern structure')
% load('roi_name.mat')

OC_pattern=struct;

for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);disp(sbj_n)
    load(fullfile('G:\JSR\240922_new_fmri\data_fmri_seg',string(strcat(string(sbj_n), '.mat'))))
    
    % set roi to number
    for hpc_or_ctx={'hpc','ctx'}
        roi_nn=seg.seg_fit.(hpc_or_ctx{:}).roi_list_nn;
        
        % set roi info
        if strcmp(hpc_or_ctx,'hpc')
            curr_mask=roi_nn(hpc_select);
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_mask=roi_nn(ctx_select);
            curr_names=roi_ctx_name(ctx_select);
        end
        sbj_nii=dir(fullfile(singletrial_path,sprintf('sub_%d',sbj_n),'*.nii'));
        
        for r=1:numel(curr_names)
            for t=1:height(sbj_nii) % 32 trials for main / 15 trials for ODT
                curr_beta=niftiread(fullfile(sbj_nii(1).folder,sbj_nii(t).name));
                pattern_roi=curr_beta(curr_mask{r});
                
                OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{1,t}=pattern_roi;
                OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,t}=obj_num{1,sbj_i}(t);%object정보
                OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{4,t}=context{1,sbj_i}{t}; %context정보
                
            end
            % if ~strcmp('main', 'main'); OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,:}=reg_1111{1, sbj_i}.ODT.(pp).object(1:12)'; end
            OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{2,r}=curr_names{r};
            OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)=roc_cell{sbj_n}';
            
        end
        OC_pattern.(hpc_or_ctx{:}).raw_roi=roi_nn;
    end
    % subject(27) - roi region(30) - trial(32)
end

save(string(fullfile(root_path,'pattern_similarity',m_o,'OC_pattern.mat')),"OC_pattern");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  REWARDING  %%%%%%%%%%%%%%%%%%%
%% rewarding!!!!
rewarding = struct;

% input
cd('G:\JSR\241104_new_fmri\pattern_similarity\main_1111')
mkdir('./rewarding');


for h_or_b={'half','block'}
    
    if strcmp(h_or_b{:},'half')
        trials=1:2;
    elseif strcmp(h_or_b{:},'block')
        trials=1:4;
    end
    
    
    for tr=trials
        for ct={'C','F'}
            contexts=ct{:};
            
            %%
            load(string(fullfile(root_path, 'pattern_similarity',m_o,'OC_pattern.mat')));
            
            for hpc_or_ctx={'hpc','ctx'}
                
                if strcmp(hpc_or_ctx,'hpc')
                    curr_names=roi_hpc_name(hpc_select);
                else
                    curr_names=roi_ctx_name(ctx_select);
                end
                bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
                
                
                for sbj_i=1:numel(sbj_id_list)
                    sbj_n = sbj_id_list(sbj_i);
                    for r=1:numel(curr_names)
                        
                        curr_s_r=OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r};
                        bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');
                        bulk_corr_ori=bulk_corr;
                        
                        if strcmp(contexts, 'F')
                            context_idx = find(cell2mat(OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(4, :))==1);
                        else
                            context_idx = find(cell2mat(OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(4, :))==2);
                        end
                        
                        corr_rej_idx = intersect(context_idx, find(cellfun(@(x) strcmp(x, 'corr_rej'), OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(2, :))));
                        hit_idx = intersect(context_idx, find(cellfun(@(x) strcmp(x, 'hit'), OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(2, :))));
                        false_idx = intersect(context_idx, find(cellfun(@(x) strcmp(x, 'false'), OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(2, :))));
                        miss_idx = intersect(context_idx, find(cellfun(@(x) strcmp(x, 'miss'), OC_pattern.(hpc_or_ctx{:}).data{1, sbj_i}{1, r}(2, :))));
                        
                        %% 진짜 corr_rej_idx만 할건지(correct only), false positive도 섞어서 할건지(corr+incorr).
                        % hit+miss / correct rejection + false positive (corr+incorr)로 할거면 아래 코드 실행하기
                        corr_rej_idx=[corr_rej_idx false_idx];
                        hit_idx=[hit_idx miss_idx];
                        
                        
                        %% all? / half? / block?
                        
                        % 1. all trials
                        % put nothing
                        
                        % 2. half 1,2
                        if strcmp(h_or_b{:},'half')
                            if tr == 1
                                corr_rej_idx = corr_rej_idx(corr_rej_idx < 17);
                                hit_idx = hit_idx(hit_idx < 17);
                            elseif tr == 2
                                corr_rej_idx = corr_rej_idx(corr_rej_idx > 16);
                                hit_idx = hit_idx(hit_idx > 16);
                            end
                            
                            % 3. block 1,2
                        elseif strcmp(h_or_b{:},'block')
                            if tr == 1
                                corr_rej_idx = corr_rej_idx(corr_rej_idx < 9);
                                hit_idx = hit_idx(hit_idx <9);
                            elseif tr == 2
                                corr_rej_idx = corr_rej_idx(corr_rej_idx > 8 & corr_rej_idx < 17);
                                hit_idx = hit_idx(hit_idx > 8 & hit_idx < 17);
                            elseif tr == 3
                                corr_rej_idx = corr_rej_idx(corr_rej_idx > 16  & corr_rej_idx < 25);
                                hit_idx = hit_idx(hit_idx > 16 & hit_idx < 25);
                            elseif tr == 4
                                corr_rej_idx = corr_rej_idx(corr_rej_idx > 24);
                                hit_idx = hit_idx(hit_idx > 24);
                            end
                        end
                        
                        
                        rewarding.hit.(contexts).(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(hit_idx,hit_idx),1);
                        rewarding.corr_rej.(contexts).(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(corr_rej_idx,corr_rej_idx),1);
                        
                        rewarding.miss.(contexts).(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(miss_idx,miss_idx),1);
                        rewarding.false.(contexts).(hpc_or_ctx{:}){r,sbj_i}=triu(bulk_corr(false_idx,false_idx),1);
                        
                        
                        
                        %% col1: hit / col2: corr_rej / col3: corr_rej-hit // row: subjects
                        rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}){sbj_i,1}=nanmean(nonzeros(rewarding.hit.(contexts).(hpc_or_ctx{:}){r,sbj_i}));
                        rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}){sbj_i,2}=nanmean(nonzeros(rewarding.corr_rej.(contexts).(hpc_or_ctx{:}){r,sbj_i}));
                        
                        rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}){sbj_i,3}=rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}){sbj_i,2}-rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}){sbj_i,1};
                        
                        rewarding.(contexts).(hpc_or_ctx{:}).idx.corr_rej{sbj_i}=corr_rej_idx;
                        rewarding.(contexts).(hpc_or_ctx{:}).idx.hit{sbj_i}=hit_idx;
                        rewarding.(contexts).(hpc_or_ctx{:}).idx.false{sbj_i}=false_idx;
                        rewarding.(contexts).(hpc_or_ctx{:}).idx.miss{sbj_i}=miss_idx;
                        
                    end
                end
                % bilateral
                ori_name=strrep(bi_names,'Bi_','');
                for r=1:numel(bi_names)
                    paat.(contexts).L=cell2mat(rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r}));
                    paat.(contexts).R=cell2mat(rewarding.(contexts).(hpc_or_ctx{:}).(curr_names{r+(numel(curr_names)/2)}));
                    
                    paat.(contexts).Bi=(paat.(contexts).L+paat.(contexts).R)/2;
                    rewarding.(contexts).(hpc_or_ctx{:}).(bi_names{r})=paat.(contexts).Bi;
                    
                    
                    % ROC로 나눠서 넣기
                    for l=1:3
                        lateral={'L_','R_','Bi_'};
                        %                rewarding.contrast.EARLY.(hpc_or_ctx{:})=[curr_names, bi_names];
                        %                rewarding.contrast.LATE.(hpc_or_ctx{:})=[curr_names, bi_names];
                        
                        pat={paat.(contexts).L,paat.(contexts).R,paat.(contexts).Bi};
                        rewarding.(contexts).ROC_pattern.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,1);
                        rewarding.(contexts).ROC_pattern.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,2);
                        rewarding.(contexts).ROC_pattern.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(:,3);
                        
                        
                        % performer group
                        
                        rewarding.(contexts).ROC_pattern.EARLY.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,1);
                        rewarding.(contexts).ROC_pattern.EARLY.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,2);
                        rewarding.(contexts).ROC_pattern.EARLY.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,3);
                        
                        rewarding.(contexts).ROC_pattern.LATE.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,1);
                        rewarding.(contexts).ROC_pattern.LATE.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,2);
                        rewarding.(contexts).ROC_pattern.LATE.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,3);
                        
                        rewarding.(contexts).ROC_pattern.FAIL.hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,1);
                        rewarding.(contexts).ROC_pattern.FAIL.corr_rej.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,2);
                        rewarding.(contexts).ROC_pattern.FAIL.subtract_corr_hit.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,3);
                        
                        
                        rewarding.(contexts).EARLY.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_early,1:2); %col 1: hit, col 2: corr_rej
                        rewarding.(contexts).LATE.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_late,1:2); %col 1: hit, col 2: corr_rej
                        rewarding.(contexts).FAIL.(hpc_or_ctx{:}).(strcat(lateral{l},ori_name{r}))=pat{l}(idx_fail,1:2); %col 1: hit, col 2: corr_rej
                        
                        rewarding.(contexts).contrast.EARLY.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_early,3))];
                        rewarding.(contexts).contrast.LATE.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_late,3))];
                        rewarding.(contexts).contrast.FAIL.(hpc_or_ctx{:})=[{strcat(lateral{l},ori_name{r})};num2cell(pat{l}(idx_fail,3))];
                        
                    end
                end
                %%
                fn = fieldnames(rewarding.(contexts).LATE.(hpc_or_ctx{:}));
                for fn_i = 1:length(fn)
                    % EARLY contrast
                    rewarding.(contexts).contrast.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).contrast.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.EARLY.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % LATE contrast
                    rewarding.(contexts).contrast.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).contrast.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.LATE.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % FAIL contrast
                    rewarding.(contexts).contrast.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).contrast.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.FAIL.subtract_corr_hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % EARLY HIT
                    rewarding.(contexts).HIT.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).HIT.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.EARLY.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % LATE HIT
                    rewarding.(contexts).HIT.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).HIT.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.LATE.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % FAIL HIT
                    rewarding.(contexts).HIT.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).HIT.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.FAIL.hit.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % EARLY CORR_REJ
                    rewarding.(contexts).CORR_REJ.EARLY.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).CORR_REJ.EARLY.(hpc_or_ctx{:})(2:length(idx_early)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.EARLY.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % LATE CORR_REJ
                    rewarding.(contexts).CORR_REJ.LATE.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).CORR_REJ.LATE.(hpc_or_ctx{:})(2:length(idx_late)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.LATE.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                    % LATE CORR_REJ
                    rewarding.(contexts).CORR_REJ.FAIL.(hpc_or_ctx{:}){1, fn_i} = fn{fn_i};
                    rewarding.(contexts).CORR_REJ.FAIL.(hpc_or_ctx{:})(2:length(idx_fail)+1, fn_i) = cellfun(@(x) x, num2cell(rewarding.(contexts).ROC_pattern.FAIL.corr_rej.(hpc_or_ctx{:}).(fn{fn_i})), 'UniformOutput', false);
                    
                end
                
            end
            
            save(strcat('./rewarding/', sprintf('%s_%d',h_or_b{:},tr) ,'_rewarding_patterns.mat'),"rewarding")
        end
        
        %% excel
        for perf={'EARLY','LATE','FAIL'}
            for reg={'hpc','ctx'}
                fn = fieldnames(rewarding.(contexts).LATE.(reg{:}));
                for hcc={'HIT','CORR_REJ','contrast'}
                    
                    T_F = cell2mat(rewarding.F.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    T_C = cell2mat(rewarding.C.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    
                    %             %both context
                    T=array2table((T_F + T_C) / 2,'VariableNames', fn');
                    writetable(T, strcat('./rewarding/', sprintf('%s_%d',h_or_b{:},tr),'_rewarding_pattern.xlsx'),'Sheet',strcat(perf{:},'_',hcc{:},'_',reg{:}));
                    
                    %             % each context
                    T=array2table(T_F,'VariableNames', fn');
                    writetable(T, strcat('./rewarding/', sprintf('%s_%d',h_or_b{:},tr),'_Forest_rewarding_pattern.xlsx'),'Sheet',strcat(perf{:},'_',hcc{:},'_',reg{:}));
                    
                    T=array2table(T_C,'VariableNames', fn');
                    writetable(T, strcat('./rewarding/', sprintf('%s_%d',h_or_b{:},tr),'_City_rewarding_pattern.xlsx'),'Sheet',strcat(perf{:},'_',hcc{:},'_',reg{:}));
                    
                end
            end
        end
    end
    
    %% half 2-1
    if strcmp(h_or_b{:},'half')
        
        half_1=load('./rewarding/half_1_rewarding_patterns.mat');
        half_2=load('./rewarding/half_2_rewarding_patterns.mat');
        
        for perf={'EARLY','LATE','FAIL'}
            for reg={'hpc','ctx'}
                fn = fieldnames(rewarding.(contexts).FAIL.(reg{:}));
                
                for hcc={'HIT','CORR_REJ','contrast'}
                    hf1_f=cell2mat(half_1.rewarding.F.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    hf1_c=cell2mat(half_1.rewarding.C.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    
                    hf2_f=cell2mat(half_2.rewarding.F.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    hf2_c=cell2mat(half_2.rewarding.F.(hcc{:}).(perf{:}).(reg{:})(2:end,:));
                    hf1=(hf1_f+hf1_c)/2;
                    hf2=(hf2_f+hf2_c)/2;
                    
                    T_cont = array2table(hf2-hf1,'VariableNames', fn');
                    writetable(T_cont, './rewarding/HALF_rewarding_pattern.xlsx','Sheet',strcat('half2-1_', perf{:},'_',hcc{:},'_',reg{:}));
                    
                    T_half1 = array2table(hf1,'VariableNames', fn');
                    writetable(T_half1, './rewarding/HALF_rewarding_pattern.xlsx','Sheet',strcat('HALF1_', perf{:},'_',hcc{:},'_',reg{:}));
                    
                    
                    T_half2 = array2table(hf2,'VariableNames', fn');
                    writetable(T, './rewarding/HALF_rewarding_pattern.xlsx','Sheet',strcat('HALF2_', perf{:},'_',hcc{:},'_',reg{:}));
                    
                end
            end
        end
        
    end
end

% rewarding structure 만들기 (총정리)
reward=struct;
for perf = {'EARLY', 'LATE', 'FAIL'}
    for hcc = {'HIT', 'CORR_REJ', 'contrast'}
        for reg = {'hpc', 'ctx'}
            
            half_1 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'half', 1), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            half_2 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'half', 2), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            block_1 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'block', 1), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            block_2 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'block', 2), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            block_3 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'block', 3), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            block_4 = readtable(strcat('./rewarding/', sprintf('%s_%d', 'block', 4), '_rewarding_pattern.xlsx'), 'Sheet', strcat(perf{:}, '_', hcc{:}, '_', reg{:}));
            
            for nm = half_1.Properties.VariableNames
                reward.half.(perf{:}).(hcc{:}).(reg{:}).(nm{:}) = [half_1.(nm{:}) half_2.(nm{:})];
                reward.block.(perf{:}).(hcc{:}).(reg{:}).(nm{:}) = [block_1.(nm{:}) block_2.(nm{:}) block_3.(nm{:}) block_4.(nm{:})];
            end
        end
    end
end
save('./rewarding/reward.mat','reward');

% rewarding Plotting
disp(pwd); mkdir('./rewarding/figures/')

for h_or_b = {'half', 'block'}
    for perf = {'EARLY', 'LATE', 'FAIL'}
        for hcc = {'HIT', 'CORR_REJ', 'contrast'}
            for reg = {'hpc', 'ctx'}
                
                fields = fieldnames(reward.(h_or_b{:}).(perf{:}).(hcc{:}).(reg{:}));
                num_fields = length(fields);
                num_plots = ceil(num_fields / 18);
                
                for plot_idx = 1:num_plots
                    figure('Position',[-1075,-105,1073,1707]);
                    sgtitle_name=sprintf('%s.%s.%s.%s', h_or_b{:}, perf{:}, strrep(hcc{:},'_','.'), reg{:});
                    sgtitle(sgtitle_name, 'FontSize', 14, 'FontWeight', 'bold');
                    
                    for i = (plot_idx - 1) * 18 + 1 : min(plot_idx * 18, num_fields)
                        subplot(6,3, i - (plot_idx - 1) * 18);
                        nm = fields{i};
                        data = reward.(h_or_b{:}).(perf{:}).(hcc{:}).(reg{:}).(nm);
                        
                        % Bar plot with error bars
                        x = 1:size(data, 2);
                        
                        means = mean(data);
                        sems = std(data) / sqrt(size(data, 1));
                        
                        bar_handle = bar(x, means, 'FaceColor', 'flat');
                        hold on;
                        
                        errorbar(x, means, sems, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
                        
                        title(strrep(nm,'_','.'), 'FontSize', 10, 'FontWeight', 'bold');
                        
                        % Adding stats
                        for j = 1:length(x)
                            [~, p] = ttest(data(:, j));
                            if p < 0.05
                                if means(j) + sems(j) >0
                                    text(x(j), means(j) + sems(j) + 0.005, '*', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                                elseif means(j) + sems(j) <=0
                                    text(x(j), 0.01, '*', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                                end
                            end
                        end
                        curr_y = ylim;y1 = curr_y(1);y2 = max(means) + max(sems) + 0.01;ylim([y1, y2]);
                        
                        box off; hold off;
                    end
                    
                    saveas(gcf, strcat('./rewarding/figures/', sgtitle_name, '_', num2str(plot_idx), '.png'));
                    %                     close(gcf);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% same Ps, diff PS
load('perform_group.mat')

same_ctxt=struct;diff_ctxt=struct;
cd(fullfile(root_path,'pattern_similarity','main'))
load('OC_pattern.mat');


for hpc_or_ctx={'hpc','ctx'}
    if strcmp(hpc_or_ctx,'hpc')
        curr_names=roi_hpc_name(hpc_select);
    else
        curr_names=roi_ctx_name(ctx_select);
    end
    
    
    for sbj_i=1:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);
        for r=1:numel(curr_names)
            
            curr_s_r=OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r};
            bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');
            
            
            f_idx = find(cellfun(@(x) x == 1, curr_s_r(4,1:4)));
            curr_f= cell2mat(OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,f_idx));
            curr_c= setdiff(4:7,curr_f);
            
            f1_idx=find(cell2mat(curr_s_r(3,:))==curr_f(1));
            f2_idx=find(cell2mat(curr_s_r(3,:))==curr_f(2));
            c1_idx=find(cell2mat(curr_s_r(3,:))==curr_c(1));
            c2_idx=find(cell2mat(curr_s_r(3,:))==curr_c(2));
            
              forest_idx=find(cellfun(@(x) x == 1, curr_s_r(4,:)));
              city_idx=find(cellfun(@(x) x == 2, curr_s_r(4,:)));
                        
              if ~isempty(intersect(forest_idx,city_idx))
            warning('%s: wrong context indexing!!!',sbj_n)
              end
            %% same
            
            bulk_same_f=bulk_corr(f1_idx,f2_idx);
            bulk_same_c=bulk_corr(c1_idx,c2_idx);
            same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[bulk_same_f;bulk_same_c];
            
            m=[];
            for lap=1:8
                corr_f=bulk_same_f;
                corr_c=bulk_same_c;
                
                same_f=corr_f(lap,lap);
                same_c=corr_c(lap,lap);
                
                m(lap)=mean([same_f;same_c],"all");
                same_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{lap,r}=m(lap);
            end
            
            for block=1:4
                same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
            end
            
            %% diff
            % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
            
            diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
            diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
            diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}=bulk_corr(c1_idx,f2_idx);
            diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);
            
            
            m=[];
            for lap=1:8
                corr_f1=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1};
                corr_f2=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2};
                corr_c1=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3};
                corr_c2=diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4};
                
                diff_f=[corr_f1(lap,lap);corr_f2(lap,lap)];
                diff_c=[corr_c1(lap,lap);corr_c2(lap,lap)];
                
                m(lap)=mean([diff_f;diff_c],"all");
                diff_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{lap,r}=m(lap);
            end
            
            for block=1:4
                diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
            end
        end
        
        %% all_mean and half
        same_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
        same_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
        diff_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
        diff_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
        
        same_ctxt.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(same_ctxt.(hpc_or_ctx{:}).half(1,sbj_i)));
        diff_ctxt.(hpc_or_ctx{:}).all_mean{1,sbj_i}= mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).half(1,sbj_i)));
    end
end
save('main_patterns.mat',"same_ctxt","diff_ctxt","sbj_id_list")


%%
now=struct;
load('main_patterns.mat');
p={'lap','block','half','all_mean'};

for diff_same={'diff','same'}
    data={};
    for hpc_or_ctx={'hpc','ctx'}
        if strcmp(hpc_or_ctx,'hpc')
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_names=roi_ctx_name(ctx_select);
        end
        
        for pp=p
            curr_pp={};
            for r=1:numel(curr_names)
                for sbj_i=1:numel(sbj_id_list)
                    temp_diff=diff_ctxt.(hpc_or_ctx{:}).(pp{:}){1,sbj_i}(:,r)';
                    temp_same=same_ctxt.(hpc_or_ctx{:}).(pp{:}){1,sbj_i}(:,r)';
                    for ts=1:numel(temp_same)
                        now.(hpc_or_ctx{:}).(pp{:}).same.(curr_names{r}){sbj_i,ts}=temp_same(ts);
                        now.(hpc_or_ctx{:}).(pp{:}).diff.(curr_names{r}){sbj_i,ts}=temp_diff(ts);
                        curr_pp{ts}=sprintf('%s_%d',pp{:},ts);
                    end
                    
                    data = now.(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(curr_names{r});
                end
                
                T = cell2table(data, 'VariableNames', curr_pp);
                writetable(T, string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet',curr_names{r});
            end
        end
    end
end

%
%% bilateral region
% cd('G:\JSR\241104_new_fmri\pattern_similarity\main')
pat=struct;

% % % 비교용: 0922
% path=('G:\JSR\241104_new_fmri\pattern_similarity\main_0922')
% cd('G:\JSR\240922_new_fmri\pattern_similarity\main\test_good');

for diff_same={'diff','same'}
    data = [];
    for hpc_or_ctx={'hpc','ctx'}
        if strcmp(hpc_or_ctx,'hpc')
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_names=roi_ctx_name(ctx_select);
        end
        bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
        
        p={'lap','block','half','all_mean'};
        for pp=p
            T_Bi=[];
            
            for r=1:numel(bi_names)
                
                T_L=readtable(string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet', curr_names{r});
                T_R=readtable(string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet', curr_names{r+(numel(curr_names)/2)});
                
% % 비교용: 0922
% T_L=readtable(string(strcat('main_within_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet', curr_names{r});
% T_R=readtable(string(strcat('main_within_',(diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_final_pattern_0922.xlsx')),'Sheet', curr_names{r+(numel(curr_names)/2)});
                            
                pat.(pp{:}).(hpc_or_ctx{:}).(diff_same{:}).(curr_names{r})=T_L;
                pat.(pp{:}).(hpc_or_ctx{:}).(diff_same{:}).(curr_names{r+(numel(curr_names)/2)})=T_R;
                
                for col=1:numel(T_L(1,:))
                    T_Bi(:,col)=mean([table2array(T_L(:,col)),table2array(T_R(:,col))],2);
                end
                
                now.(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(bi_names{r})=T_Bi;
                
                T = array2table(T_Bi, 'VariableNames',T_L.Properties.VariableNames);
                writetable(T, string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet',bi_names{r});
                pat.(pp{:}).(hpc_or_ctx{:}).(diff_same{:}).(bi_names{r})=T;
                
                % performance group
                
                % fail group
                [pat]=func_save_pat(T, T_L, T_R, idx_fail, 'FAIL', diff_same{:}, hpc_or_ctx{:}, pp{:}, bi_names, curr_names, r,pat);
                % good learner
                idx_good = sort([idx_early idx_late]);
                [pat]=func_save_pat(T, T_L, T_R, idx_good, 'GOOD', diff_same{:}, hpc_or_ctx{:}, pp{:}, bi_names, curr_names, r,pat);
                % early learner
                [pat]=func_save_pat(T, T_L, T_R, idx_early, 'EARLY', diff_same{:}, hpc_or_ctx{:}, pp{:}, bi_names, curr_names, r,pat);
                % late learner
                [pat]=func_save_pat(T, T_L, T_R, idx_late, 'LATE', diff_same{:}, hpc_or_ctx{:}, pp{:}, bi_names, curr_names, r,pat);
            end
        end
    end
end
final_pat=pat;
save('now.mat','now')
save('final_pat.mat','final_pat','idx_good','idx_early','idx_late','idx_fail')

%% in main, same-diff
for hpc_or_ctx={'hpc','ctx'}
    if strcmp(hpc_or_ctx,'hpc')
        curr_names=roi_hpc_name(hpc_select);
    else
        curr_names=roi_ctx_name(ctx_select);
    end
    bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
    roi=[curr_names,bi_names];
    
    
    p={'lap','block','half','all_mean'};
    for pp=p
        
        temp_tbl=string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx'));
        sheets=sheetnames(temp_tbl);
        for shts=1:numel(sheets)
            same_ptn=readtable(string(strcat('same_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet',sheets(shts));
            diff_ptn=readtable(string(strcat('diff_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet',sheets(shts));
            
            s_d_ptn=table2array(same_ptn)-table2array(diff_ptn);
            
            T = array2table(s_d_ptn, 'VariableNames', same_ptn.Properties.VariableNames);
            writetable(T,string(strcat('same-diff_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet',sheets(shts));
            pat.same_diff.(hpc_or_ctx{:}).(pp{:})=T;

        end
    end
end
final_pat=pat;
save('final_pat.mat','final_pat','idx_good','idx_early','idx_late','idx_fail')

%% plot

disp(pwd); mkdir('./figures/')
func_pat_plot(final_pat.perform)


%% ROI only
load('final_pat.mat')
pat=final_pat;
perf_group = fieldnames(pat.perform)';


for pp=perf_group
    for diff_same={'diff','same'}
        for roi_name={'Bi_DGCA3','Bi_CA1'}
            T=pat.perform.(pp{:}).block.hpc.(diff_same{:}).(roi_name{:});
            writetable(T,'ROI_block_PS.xlsx','Sheet',strcat(pp{:},'_',diff_same{:},'_',roi_name{:}))
        end
    end
end
            



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
function func_pat_plot(pat_in)

group_names = fieldnames(pat_in);
for g = 1:length(group_names)
    group_name = group_names{g};
    pp_names = fieldnames(pat_in.(group_name));
    
    for p = 1:length(pp_names)
        pp = pp_names{p};
        hpc_ctx_names = fieldnames(pat_in.(group_name).(pp));
        
        for h = 1:length(hpc_ctx_names)
            hpc_or_ctx = hpc_ctx_names{h};
            diff_same_names = fieldnames(pat_in.(group_name).(pp).(hpc_or_ctx));
            
            for d = 1:length(diff_same_names)
                diff_same = diff_same_names{d};
                
                fields = fieldnames(pat_in.(group_name).(pp).(hpc_or_ctx).(diff_same));
                num_fields = length(fields);
                num_plots = ceil(num_fields / 18);
                
                for plot_idx = 1:num_plots
                    figure('Position',[-1075,-105,1073,1707]);
                    sgtitle_name = sprintf('%s.%s.%s.%s',  pp, hpc_or_ctx, group_name, diff_same);
                    sgtitle(strrep(sgtitle_name,'_','.'), 'FontSize', 14, 'FontWeight', 'bold');
                    
                    for i = (plot_idx - 1) * 18 + 1 : min(plot_idx * 18, num_fields)
                        subplot(6, 3, i - (plot_idx - 1) * 18);
                        var_name = fields{i};
                        data = table2array(pat_in.(group_name).(pp).(hpc_or_ctx).(diff_same).(var_name));
                        
                        % Bar plot with error bars
                        x = 1:size(data, 2);
                        means = mean(data);
                        sems = std(data) / sqrt(size(data, 1));
                        
                        bar_handle = bar(x, means, 'FaceColor', 'flat');
                        hold on;
                        
                        errorbar(x, means, sems, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
                        curr_y = ylim; y1 = curr_y(1);y2 = max(means + sems) + 0.01; ylim([y1, y2]);
                        
                        title(strrep(var_name,'_','.'), 'FontSize', 10, 'FontWeight', 'bold');
                        
                        % Adding stats
                        for j = 1:length(x)
                            [~, p_val] = ttest(data(:, j));
                            if p_val < 0.05
                                 if means(j) + sems(j) >0
                                    text(x(j), means(j) + sems(j) + 0.005, '*', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                                elseif means(j) + sems(j) <=0
                                    text(x(j), 0.01, '*', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                                end
                            end
                        end
                        
                        hold off;
                    end
                    
                    saveas(gcf, strcat('./figures/', sgtitle_name, '_', num2str(plot_idx), '.png'));
                                        close(gcf);
                end
            end
        end
    end
end

end



%% function
function [final_pat]=func_save_pat(T, T_L, T_R, idx_group, group_name, diff_same, hpc_or_ctx, pp, bi_names, curr_names, r, pat)
writetable(T(idx_group,:), string(strcat(group_name, '_', diff_same, '_', hpc_or_ctx, '_', pp, '_pattern.xlsx')), 'Sheet', bi_names{r});
writetable(T_L(idx_group,:), string(strcat(group_name, '_', diff_same, '_', hpc_or_ctx, '_', pp, '_pattern.xlsx')), 'Sheet', curr_names{r});
writetable(T_R(idx_group,:), string(strcat(group_name, '_', diff_same, '_', hpc_or_ctx, '_', pp, '_pattern.xlsx')), 'Sheet', curr_names{r + (numel(curr_names) / 2)});
pat.perform.(group_name).(pp).(hpc_or_ctx).(diff_same).(bi_names{r})=T(idx_group,:);
pat.perform.(group_name).(pp).(hpc_or_ctx).(diff_same).(curr_names{r})=T_L(idx_group,:);
pat.perform.(group_name).(pp).(hpc_or_ctx).(diff_same).(curr_names{r + (numel(curr_names) / 2)})=T_R(idx_group,:);
final_pat=pat;
end
