clc;clear
%% single trial naming
% zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
% addpath(zcode_path)
% root_path = 'C:\Users\User\Desktop\JSR';%90번컴
root_path = 'G:\JSR\240922_new_fmri'; %69번컴
cd(fullfile(root_path, 'spm_prep_glm_0922')); addpath('G:\JSR\sbj39to41')

load('../bhv/data_bhv_log_table/total/num_sbj_events.mat')
load('../regressors_GLM_0922.mat')
load('../bhv/learn_curv_info.mat');
roi_hpc_name=strrep(roi_hpc_name,'CA23DG','DGCA3');


% sbj_id_list_good=setdiff(sbj_id_list_41,[2,3,5,6,7,15,21,22,30,32,33,38]);
sbj_id_list=sbj_id_list_41;
tbl=num_sbj_events;
%% setting
main_or_ODT = {'pre-ODT','post-ODT'}; % {'main','pre-ODT','post-ODT'} for for loops


hpc_select=[1:5,7:17,19:24];
ctx_select=1:30;

%%
for m_o = main_or_ODT
    disp(m_o)
    
    singletrial_path = fullfile(root_path, 'spm_prep_glm_0922',m_o{:}, 'single_trial');
    
    % % 잘못만들었을때 delete할 수 있는 코드
    %     for sbj_i=1:numel(sbj_id_list)
    %                 sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
    %     files=dir(fullfile(singletrial_path,sbj_n, '*_0*'));
    %     for i=1:height(files)
    %     delete(fullfile(files(1).folder,files(i).name))
    %     end
    %     end
    %     end
    % % %
    
    %% rename_ trial
    if strcmp(m_o{:}, 'main')
        disp('rename_ trial')
        
        for sbj_i=1:numel(sbj_id_list)
            sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
            pp='obj_show';
            sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001',strcat(pp,'_*')));
            roc=cellfun(@(x) extractBetween(x,strcat(pp,'_'),'_0'),{sbj_OC.name});
            d=reg_0922{1, sbj_i}.trial_detail;
            un=unique(roc);
            
            sidx=[];
            for r=1:numel(un)
                idx=d.(un{r})';
                sidx=[sidx idx];
            end
            [b,iii]=sort(sidx);
            if max(iii)==32
                tbl.roc(tbl.session==sbj_id_list(sbj_i))=roc(iii);
                for i=1:height(sbj_OC)
                    beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
                    beta_rename=sprintf('%.2d_%s.nii',sidx(i), sbj_OC(i).name);
                    copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));
                end
            else
                break;
                error('The total number of trials is not 32!');
            end
        end
        
    else
        for sbj_i=1:numel(sbj_id_list)
            sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));disp(sbj_n)
            pp=extractBefore(m_o{:},'-');
            sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001',strcat(pp,'*')));
            obj=cellfun(@(x) extractBetween(x,strcat(pp,'_'),'_0'),{sbj_OC.name},'UniformOutput',false);
            %         obj=repelem(["f1","f2","c1","c2","target"],3);
            temp_onset=reg_0922{1, sbj_i}.ODT.(pp).regress_onset;
            if contains(obj{1},'city')
                temp_onset=[temp_onset(3:4); temp_onset(1:2); temp_onset(5)]; %["f1","f2","c1","c2","target"] regressor 순서로 바꿔줌
            end
            [tr,tr_idx]=sort(cell2mat(temp_onset));
            
            
            beta_rename={};
            for tridx=1:length(tr_idx)
                beta_rename{tridx}=sprintf('%.2d_%s',tridx,sbj_OC(tr_idx(tridx)).name);
                beta=fullfile(sbj_OC(tr_idx(tridx)).folder,sbj_OC(tr_idx(tridx)).name,'beta_0001.nii');
                copyfile(beta, fullfile(singletrial_path,sbj_n, strcat(beta_rename{tridx},'.nii')));
            end
            
            
            reg_0922{1, sbj_i}.ODT.(pp).object=obj(tr_idx);
            reg_0922{1, sbj_i}.ODT.(pp).object_idx=tr_idx;
        end
        save('regressors_GLM_0922.mat',"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name')
    end
    
    
    
    %% Trial detail
    % load("sbj_events.mat")
    
    if strcmp(m_o{:}, 'main')
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
            
            startIdx = (sbj_id_list(sbj_i)-1)*32 + 1;
            endIdx = min(sbj_id_list(sbj_i)*32, height(tbl));
            obj_num{sbj_i} = tbl.Object(startIdx:endIdx)';
            for cont=1:32
                if tbl.Association(cont) ==1
                    ctxt{cont} = tbl.Context(cont);
                elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==1
                    ctxt{cont} = 2;
                elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==2
                    ctxt{cont} = 1;
                end
                context{sbj_i}=ctxt;
            end
            
            reg_0922{1, sbj_i}.trial_detail.context=cell2mat(context{sbj_i});
            reg_0922{1, sbj_i}.trial_detail.object=obj_num{sbj_i};
        end
        save('regressors_GLM_0922.mat',"reg_0922","sbj_id_list_41",'roi_ctx_name','roi_hpc_name')
        
    end
    
    
    %% pattern structure
    disp('pattern structure')
    % load('roi_name.mat')
    
    OC_pattern=struct;
    
    singletrial_path = fullfile(root_path, 'spm_prep_glm_0922',m_o{:}, 'single_trial');
    
    
    for sbj_i=1:numel(sbj_id_list)
        sbj_n = sbj_id_list(sbj_i);disp(sbj_n)
        load(fullfile(root_path,'data_fmri_seg',string(strcat(string(sbj_n), '.mat'))))
        %         if sbj_i==1
        %             roi_hpc_name=seg.seg_fit.hpc.roi_name_list;
        %             roi_ctx_name=seg.seg_fit.ctx.roi_name_list;
        %         end
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
                    OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{2,t}=extractBetween(sbj_nii(t).name,strcat(pp,'_'),'_0');
                    if strcmp(m_o{:}, 'main')
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,t}=obj_num{1,sbj_i}(t);%object정보
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{4,t}=context{1,sbj_i}{t}; %context정보
                    else
                        OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,t}=reg_0922{1, sbj_i}.ODT.(pp).object{t};
                    end
                end
                % if ~strcmp(m_o{:}, 'main'); OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}{3,:}=reg_0922{1, sbj_i}.ODT.(pp).object(1:12)'; end
                OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{2,r}=curr_names{r};
                
            end
            OC_pattern.(hpc_or_ctx{:}).raw_roi=roi_nn;
        end
        % subject(27) - roi region(30) - trial(32)
    end
    
    save(string(fullfile(root_path,'spm_prep_glm_0922',m_o,'OC_pattern.mat')),"OC_pattern");
    
    
end
% roi_ctx_name=strrep(roi_ctx_name,'.','_');
% roi_hpc_name=strrep(roi_hpc_name,'.','_');
% save('roi_name.mat','roi_ctx_name','roi_hpc_name')

% %%%%%%%%%%%%%%%%%%% 이 아래 OC pattern에 맞춰서 다시 수정하기!!!
% % rewarding!!!!
% main_or_ODT = {'pre-ODT','post-ODT'};
% load('regressors_GLM_0922.mat')
%
%
%     load(string(strcat('C:\Users\User\Desktop\JSR\spm_prep_glm_0922\',main_or_ODT,'\OC_pattern.mat')));
%
%
%
% if strcmp(hpc_or_ctx,'hpc')
%
%     curr_names=roi_hpc_name(hpc_select);
% else
%
%     curr_names=roi_ctx_name(ctx_select);
% end
%
%
% rewarding=struct;
% for sbj_i=1:numel(sbj_id_list)
%     sbj_n = sbj_id_list(sbj_i);
%     for r=1:numel(curr_names)
%
%         curr_s_r=OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r};
%         bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');
%         bulk_corr_ori=bulk_corr;
%
%         corr_rej_idx=find(cellfun(@(x) strcmp(x,'corr_rej'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
%         hit_idx=find(cellfun(@(x) strcmp(x,'hit'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
%         false_idx=find(cellfun(@(x) strcmp(x,'false'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
%         miss_idx=find(cellfun(@(x) strcmp(x,'miss'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(2,:)));
%
%         rewarding.hit{r,sbj_i}=triu(bulk_corr(hit_idx,hit_idx),1);
%         rewarding.corr_rej{r,sbj_i}=triu(bulk_corr(corr_rej_idx,corr_rej_idx),1);
%         rewarding.(curr_names{r}){sbj_i,1}=mean(nonzeros(rewarding.hit{r,sbj_i}));
%         rewarding.(curr_names{r}){sbj_i,2}=mean(nonzeros(rewarding.corr_rej{r,sbj_i}));
%
%         rewarding.(curr_names{r}){sbj_i,3}=rewarding.(curr_names{r}){sbj_i,2}-rewarding.(curr_names{r}){sbj_i,1};
%
%         rewarding.idx.corr_rej{sbj_i}=corr_rej_idx;
%         rewarding.idx.hit{sbj_i}=hit_idx;
%         rewarding.idx.false{sbj_i}=false_idx;
%         rewarding.idx.miss{sbj_i}=miss_idx;
%
%     end
% end
% save('rewaarding.mat',"rewarding")

% %
%% same Ps, diff PS
% load('roi_name.mat')

for m_o = main_or_ODT
    disp(m_o)
    
    same_ctxt=struct;diff_ctxt=struct;
    
    load(string(fullfile(root_path,'spm_prep_glm_0922',m_o,'OC_pattern.mat')));
    
    
    % load('PS_basic_info.mat');
    % sbj_id_list=sbj_id_list(1:28)
    
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
                
                if strcmp(m_o{:}, 'main')
                    
                    f_idx = find(cellfun(@(x) x == 1, curr_s_r(4,1:4)));
                    curr_f= cell2mat(OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,f_idx));
                    curr_c= setdiff(4:7,curr_f);
                    
                    f1_idx=find(cell2mat(curr_s_r(3,:))==curr_f(1));
                    f2_idx=find(cell2mat(curr_s_r(3,:))==curr_f(2));
                    c1_idx=find(cell2mat(curr_s_r(3,:))==curr_c(1));
                    c2_idx=find(cell2mat(curr_s_r(3,:))==curr_c(2));
                    
                else
                    f1_idx = find(cellfun(@(x) contains(x,'forest_1'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    f2_idx = find(cellfun(@(x) contains(x,'forest_2'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    c1_idx = find(cellfun(@(x) contains(x,'city_1'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    c2_idx = find(cellfun(@(x) contains(x,'city_2'),OC_pattern.(hpc_or_ctx{:}).data{1,sbj_i}{1,r}(3,:)));
                    
                end
                %% same
                %                 if issymmetric(bulk_corr)
                %                     tri_bulk=tril(bulk_corr,-1);
                %                 end
                
                %                 same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];
                beta_same=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];
                beta_same=reshape(beta_same,numel(beta_same),1);
                tri_bulk_corr=triu(bulk_corr,1);
                same_same= [tri_bulk_corr(f1_idx,f1_idx); tri_bulk_corr(c1_idx,c1_idx);tri_bulk_corr(f2_idx,f2_idx); tri_bulk_corr(c2_idx,c2_idx)];
                same_same=same_same(same_same~=0);
                same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}=[beta_same; same_same];
                
                if strcmp(m_o{:}, 'main')
                    m=[];
                    for lap=1:8
                        m(lap)=mean([same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}(lap,lap); same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}(lap+8,lap)],"all");
                    end
                    same_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m';
                    for block=1:4
                        same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
                    end
                end
                %% diff
                % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
                
                diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
                diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
                diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}=bulk_corr(c1_idx,f2_idx);
                diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);
                
                
                if strcmp(m_o{:}, 'main')
                    m=[];
                    for lap=1:8
                        m(lap)=mean([diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}(lap,lap); diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}(lap,lap);diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}(lap,lap);diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4}(lap,lap)],"all");
                    end
                    diff_ctxt.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m';
                    for block=1:4
                        diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
                    end
                end
            end
            if strcmp(m_o{:}, 'main')
                same_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                same_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                diff_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                diff_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
            end
        end
        
        save(string(strcat(m_o{:},'_patterns_0922')),"same_ctxt","diff_ctxt","sbj_id_list_41")
        
    end
end


%%

for m_o = main_or_ODT
    % main_or_ODT = {'main'};
    disp(m_o{:});
    
    
    load(string(strcat(m_o,'_patterns_0922')));
    
    for hpc_or_ctx={'hpc','ctx'}
        if strcmp(hpc_or_ctx,'hpc')
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_names=roi_ctx_name(ctx_select);
        end
        for r=1:numel(curr_names)
            for sbj_i=1:numel(sbj_id_list)
                m2=[];
                % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2
                m2=mean([mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{1}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{2}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{3}); mean(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}{4})]);
                diff_ctxt.new_mean.(hpc_or_ctx{:}).lap{1,sbj_i}{:,r}=m2;
                
                if strcmp(m_o{:}, 'main')
                    for block=1:4
                        diff_ctxt.new_mean.(hpc_or_ctx{:}).block{1,sbj_i}{block,r}=mean(m2(block*2-1:block*2));
                    end
                    diff_ctxt.new_mean.(hpc_or_ctx{:}).half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(1:2,:)));
                    diff_ctxt.new_mean.(hpc_or_ctx{:}).half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).block{1,sbj_i}(3:4,:)));
                end
            end
        end
    end
    
    save(string(strcat(m_o,'_patterns_0922')),"same_ctxt","diff_ctxt","sbj_id_list_41")
    
end

%%

now=struct;
sbj_id_list=sbj_id_list_41;

for m_o = main_or_ODT
    
    
    load(string(strcat(m_o,'_patterns_0922')));
    
    
    if strcmp(m_o{:}, 'main')
        pp='main';
    else
        pp=extractBefore(m_o{:},'-');
    end
    
    for diff_same={'diff','same'}
        data = [];
        for hpc_or_ctx={'hpc','ctx'}
            if strcmp(hpc_or_ctx,'hpc')
                curr_names=roi_hpc_name(hpc_select);
            else
                curr_names=roi_ctx_name(ctx_select);
            end
            
            for sbj_i=1:numel(sbj_id_list)
                sbj_n = sbj_id_list(sbj_i);
                for r=1:numel(curr_names)
                    if strcmp(m_o{:}, 'main')
                        now.(hpc_or_ctx{:}).diff.(curr_names{r}){sbj_i,:}=diff_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(:,r)';
                        now.(hpc_or_ctx{:}).same.(curr_names{r}){sbj_i,:}=same_ctxt.(hpc_or_ctx{:}).half{1,sbj_i}(:,r)';
                        
                        % now.diff.(curr_names{r}){sbj_i,:}=cell2mat(diff_ctxt.block  {1,sbj_i}(:,r)');
                        % now.same.(curr_names{r}){sbj_i,:}=cell2mat(same_ctxt.block{1,sbj_i}(:,r)');
                    else
                        now.(hpc_or_ctx{:}).(pp).diff.(curr_names{r}){sbj_i,:}=mean(cell2mat(diff_ctxt.(hpc_or_ctx{:}).all{r,sbj_i}),"all");
                        now.(hpc_or_ctx{:}).(pp).same.(curr_names{r}){sbj_i,:}=mean(same_ctxt.(hpc_or_ctx{:}).all{r,sbj_i},"all");
                    end
                    
                    data(sbj_i, r) = now.(hpc_or_ctx{:}).(pp).(diff_same{:}).(curr_names{r}){sbj_i};
                end
            end
            T = array2table(data, 'VariableNames', curr_names);
            if strcmp(m_o{:},'main')
                writetable(T, string(strcat(root_path, '\pattern_similarity\main\', m_o,'_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
            else
                writetable(T, string(strcat(root_path, '\pattern_similarity\ODT\', m_o,'_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
            end
        end
    end
end

%% bilateral region

for m_o = main_or_ODT
    
    for diff_same={'diff','same'}
        data = [];
        for hpc_or_ctx={'hpc','ctx'}
            if strcmp(hpc_or_ctx,'hpc')
                curr_names=roi_hpc_name(hpc_select(1:end/2));
            else
                curr_names=roi_ctx_name(ctx_select(1:end/2));
            end
            curr_names=strrep(curr_names,'Lt','Bi');
            
            
            T=readtable(string(strcat(root_path, '\pattern_similarity\ODT\', m_o,'_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
            
            
            for r=1:numel(curr_names)
                data(:, r) = mean([table2array(T(:,r)),table2array(T(:,r+numel(curr_names)))],2);
            end
            
            T_bi = array2table(data, 'VariableNames', curr_names);
            T_combi = [T, T_bi];
            
            if strcmp(m_o{:},'main')
                writetable(T_combi, string(strcat(root_path, '\pattern_similarity\main\', m_o,'_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
            else
                writetable(T_combi, string(strcat(root_path, '\pattern_similarity\ODT\', m_o,'_pattern_0922.xlsx')),'Sheet', strcat((diff_same{:}),'_',(hpc_or_ctx{:})));
            end
        end
    end
    
end
%% if ODT, post-pre
if ~strcmp(m_o{:},'main')
    
    for hpc_or_ctx={'hpc','ctx'}
        if strcmp(hpc_or_ctx,'hpc')
            curr_names=roi_hpc_name(hpc_select);
        else
            curr_names=roi_ctx_name(ctx_select);
        end
        bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
        roi=[curr_names,bi_names];
        
        for diff_same={'diff','same'}
            
            odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).subt=odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).post - odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).pre;
            
            T = array2table(odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).subt, 'VariableNames', roi);
            
            writetable(T, [root_path '\pattern_similarity\ODT\post-pre_ODT_pattern_0922.xlsx'], 'Sheet', strcat((diff_same{:}), '_', (hpc_or_ctx{:})));
            
        end
        PSdiff=odt_ptn.(hpc_or_ctx{:}).same.subt - odt_ptn.(hpc_or_ctx{:}).diff.subt;
        T = array2table(PSdiff, 'VariableNames', roi);
        
        writetable(T, [root_path '\pattern_similarity\ODT\post-pre_ODT_pattern_0922.xlsx'], 'Sheet', strcat('same-diff_', (hpc_or_ctx{:})));
    end
end

%% early/late/fail grouping
if contains(main_or_ODT,'ODT')
    if contains(main_or_ODT,'main')
        m_o_O = {'post-ODT','pre-ODT','post-pre_ODT', 'main'}; % 'main', 'post-ODT','pre-ODT','post-pre_ODT'
    else
        m_o_O = {'post-ODT','pre-ODT','post-pre_ODT'};
    end
else
    m_o_O = {'main'};
end

for m_o = m_o_O 
    if strcmp(m_o{:},'main')
        cd([root_path '\pattern_similarity\main'])
    else
        cd([root_path '\pattern_similarity\ODT'])
    end
    curr_excel=strcat(m_o{:}, '_pattern_0922.xlsx');
    shts = sheetnames(curr_excel);
    
    for sh=1:numel(shts)
        curr_sheet=shts(sh);
        T=readtable(curr_excel,'Sheet', shts(sh));
        
        T_early=T(idx_early,:);
        T_late=T(idx_late,:);
        T_fail=T(idx_fail,:);
        
        writetable(T_early, strcat('./perform_group/early_',curr_excel),'Sheet', shts(sh));
        writetable(T_late, strcat('./perform_group/late_',curr_excel),'Sheet', shts(sh));
        writetable(T_fail, strcat('./perform_group/fail_',curr_excel),'Sheet', shts(sh));
        
    end
end



% %% plot!!
%
% for m_o = main_or_ODT
%
%     for diff_same={'diff','same'}
%         data = [];
%         for hpc_or_ctx={'hpc','ctx'}
%
%                 T=readtable(string(strcat(m_o,'_pattern_0922.xlsx')),'Sheet', strcat(diff_same{:},'_',hpc_or_ctx{:}));
%
%             if strcmp(hpc_or_ctx,'hpc')
%                 curr_names=roi_hpc_name(hpc_select);
%             else
%                 curr_names=roi_ctx_name(ctx_select);
%             end
%             bi_names=strrep(curr_names(1:end/2),'Lt','Bi');
%
%             roi=[curr_names,bi_names];
%             roi=strrep(roi,'_','-');
%             data = table2array(T);
%
%             % if ODT, post-pre
%             if strcmp(m_o{:}, 'main')
%                 break;
%             else
%                 pp=extractBefore(m_o{:},'-');
%                 odt_ptn.(hpc_or_ctx{:}).(diff_same{:}).(pp)=data;
%             end
%
%
%             % Plotting
%
%             figure('position',[472,604,2063,894]);
%             hold on
%             bar(categorical(roi),mean(data));
%             title([m_o{:}, ': ',diff_same{:},'(',hpc_or_ctx{:},')'],'fontsize',20);
%             xlabel('Regions','fontsize',20);
%             ylabel('pattern similarity','fontsize',20);
%             ylim([0 0.3])
%             hold off
%
%                 plot_path=fullfile(root_path,'spm_prep_glm_0922','final2',m_o{:},'plot');
%
%             mkdir(plot_path)
%             saveas(gcf,[plot_path '\' diff_same{:} '_' hpc_or_ctx{:} '.png'])
%         end
%     end
% end








% %
% % % --> 이거로 69번 컴에서 graph pad prism 으로 plot그리기 (block별, half 별)
% %%
