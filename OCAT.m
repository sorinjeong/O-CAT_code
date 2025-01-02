close all; clc; clear
path_root='G:\JSR\241104_new_fmri\activity\princeton_searchlight';
addpath(genpath(path_root))
cd([path_root '\princeton-mvpa-toolbox-master'])
addpath('C:\Users\leelab\Documents\MATLAB\spm12');

mvpa_add_paths;
help tutorial_easy_spm
searchlight=struct;good_data=[];

%'./core/preproc/'에 있음
mex './core/preproc/compute_xcorr.c'

%% load OCAT data
addpath(genpath('G:\JSR\241104_new_fmri\glm_new_1111\main\single_trial'))
load('G:\JSR\241104_new_fmri\regressor_1111.mat')

% % learning timing (lap, block)
% load('G:\JSR\241104_new_fmri\learn_curv_info.mat')
% learn=struct;learn.trial=learn_trial;
% learn.lap = ceil(learn.trial/ 4);
% learn.block= ceil(learn.lap/ 2);
% save('G:\JSR\241104_new_fmri\learn_timing.mat','learn');

load('G:\JSR\241104_new_fmri\learn_timing.mat')
%         load('G:\JSR\241104_new_fmri\learn_curv_info.mat')

%% input!!
version= 'ver1'; % ver 1(forest,city) / 2(obj ID) / 3(ROC)
tm= 'lap'; % run : lap, block, trial


%% import OCAT data
for v=1:3
    version=sprintf('ver%d',v);
        for t=1:3
%     t=3;
    tt={'lap','block','trial'};
    tm=tt{t};
    
    searchlight.(version).(tm).good_data=[];searchlight.wo_anova.(version).(tm).good_data=[];
    for sbj_i=1:numel(sbj_id_list)
        % sbj_i=1
        % searchlight version 1
        % associated context
        condname.ver1={'Asso_F','Asso_C'};
        reg.ver1=double([(reg_1111{1, sbj_i}.trial_detail.context==1);  (reg_1111{1, sbj_i}.trial_detail.context==2)]);
        
        % searchlight version 2
        % object id
        condname.ver2={'Obj_1','Obj_2','Obj_3','Obj_4'};
        reg.ver2=double([(reg_1111{1, sbj_i}.trial_detail.object==4); (reg_1111{1, sbj_i}.trial_detail.object==5); ...
            (reg_1111{1, sbj_i}.trial_detail.object==6); (reg_1111{1, sbj_i}.trial_detail.object==7)]);
        
        % searchlight version 3
        % hit/miss/correct rejection / false positive
        condname.ver3={'Roc_hit','Roc_miss','Roc_corr','Roc_fals'};
        reg.ver3=double([strcmp(reg_1111{1, sbj_i}.trial_detail.all.ch,'hit'); strcmp(reg_1111{1, sbj_i}.trial_detail.all.ch,'miss'); ...
            strcmp(reg_1111{1, sbj_i}.trial_detail.all.ch,'corr_rej'); strcmp(reg_1111{1, sbj_i}.trial_detail.all.ch,'false')]);
        
        % set timepoint(run)
        timepoint.lap = repelem(1:8, 4);
        timepoint.block = repelem(1:4, 8);
        timepoint.trial = [1:32];
        
        %% 1) Initializing the subject (''subj'') structure
        subj = init_subj('ocat',sprintf('sub_%d',sbj_id_list(sbj_i)));
        summarize(subj);
        
        %% 2) Ventral temporal with GLM mask
        
        pth='G:\JSR\241104_new_fmri\activity\princeton_searchlight\OCAT_mask';
        cd(pth)
        % %double mask로부터 nifti file 생성
        % for sbj_i=1:numel(sbj_id_list)
        %
        % load(fullfile('G:\JSR\240922_new_fmri\data_fmri_seg',string(strcat(string(sbj_i), '.mat'))))
        %
        % hpcNphg = seg.seg_fit.hpc.roi_list_nn{1,25};
        % for i = 26:36
        %     hpcNphg = hpcNphg | seg.seg_fit.hpc.roi_list_nn{1,i};
        % end
        %
        % header=seg.seg_fit.header;
        % header.Datatype='double';
        % roi_mask=double(hpcNphg);
        % niftiwrite(roi_mask,fullfile(pth,sprintf('sub-%d.nii',sbj_id_list(sbj_i))),header)
        % end
        
        subj = load_spm_mask(subj,'seg_hpc_25to36',sprintf('sub-%d.nii',sbj_id_list(sbj_i)));
        summarize(subj)
        mask = get_object(subj,'mask','seg_hpc_25to36')
        
        
        
        %% 3) EPI pattern
        
        pth_lss='G:\JSR\241104_new_fmri\glm_new_1111\main\single_trial';
        pth_lss_sbj=fullfile(pth_lss, sprintf('sub_%d',sbj_id_list(sbj_i)));
        cd(pth_lss_sbj)
        
        d=dir(fullfile(pth_lss_sbj,'*.nii'))
        
        for i=1:32
            raw_filenames{i} = fullfile(d(1).folder,d(i).name);
        end
        
        subj = load_spm_pattern(subj,'epi','seg_hpc_25to36',raw_filenames);
        pats = get_object(subj,'pattern','epi')
        
        
        %% 4) Condition regressors
        curr_reg=reg.(version);
        
        subj = init_object(subj,'regressors','conds');
        whos curr_reg
        
        subj = set_mat(subj,'regressors','conds',curr_reg);
        subj = set_objfield(subj,'regressors','conds','condnames',condname.(version));
        condsobj= get_object(subj,'regressors','conds')
        
        
        %% 5) Runs selector
        subj = init_object(subj,'selector','runs');
        
        
        subj = set_mat(subj,'selector','runs', timepoint.(tm));
        runs = get_object(subj,'selector','runs')
        
        % % Try examining the contents of the ''runs'' selector matrix yourself
        % %(the 1x1210 ''mat'' field) as before, using:
        
        % runs = get_mat(subj,'selector','runs');
        
        % %You should see that it contains 1210 integer labels, one for each TR,
        % %identifying which scanning run that TR came from.
        
        summarize(subj)
        
        %% 6) Zscoring
        % open('duplicate_object.m') 한 후 datetime(true) --> datetime('now')로 변경!!
        
        subj = zscore_runs(subj,'epi','runs');
        summarize(subj,'objtype','pattern')
        
        before_zscore_mat = get_mat(subj,'pattern','epi');
        after_zscore_mat = get_mat(subj,'pattern','epi_z');
        
        
        %% 7) Creating the cross-validation indices
        subj = create_xvalid_indices(subj,'runs');
        summarize(subj,'objtype','selector')
        % % You can even peer at each selector individually, e.g.
        % xval1 = get_object(subj,'selector','runs_xval_1')
        
        
        %% 8) ANOVA
        subj = feature_select(subj,'epi_z','conds','runs_xval')
        
        summarize(subj,'display_groups',false)
        summarize(subj,'objtype','mask')
        
        
        %% 9) Cross-validation classification
        % open('cross_validation.m') 한 후 datetime(true) --> datetime('now')로 변경!!
        disp(sbj_id_list(sbj_i))
        
        class_args.train_funct_name = 'train_bp';
        class_args.test_funct_name = 'test_bp';
        class_args.nHidden = 0;
        
        [subj results] = cross_validation(subj,'epi_z','conds','runs_xval','epi_z_thresh0.05',class_args);
        if t==3
            [subj_wo_anova results_wo_anova] = cross_validation(subj,'epi_z','conds','runs_xval','seg_hpc_25to36',class_args);
        end
        
        %% 결과 확인
        %total_perf 확인 (iteration 전체의 평균값)
        results
        % acts, perfmet 확인
        % results.iterations(1)
        %guesses에서 예측한 결과를 확인할 수 있고, 그 퍼포먼스는 아래 코드로 확인 가능
        % results.iterations(1).perfmet.perf
        
        
        %% performance of pattern classification algorithms from princeton mvpa toolbox
        
        p=[results.iterations(:).perf];
        
        searchlight.(version).(tm).perf(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), p}
        if isnan(learn.(tm)(sbj_i))
            searchlight.(version).(tm).perf_learn(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), NaN}
            
        else
            %     find(p==max(p))-learn.(tm)(sbj_i)
            searchlight.(version).(tm).good_data(end+1,:)=p;
            searchlight.(version).(tm).perf_learn(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), p(learn.(tm)(sbj_i))}
        end
        searchlight.(version).(tm).perf_total(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), [results.total_perf]}
        
        % trials without ANOVA
        if t==3
            
            p=[results_wo_anova.iterations(:).perf];
            
            searchlight.wo_anova.(version).(tm).perf(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), p}
            if isnan(learn.(tm)(sbj_i))
                searchlight.wo_anova.(version).(tm).perf_learn(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), NaN}
                
            else
                %     find(p==max(p))-learn.(tm)(sbj_i)
                searchlight.wo_anova.(version).(tm).good_data(end+1,:)=p;
                searchlight.wo_anova.(version).(tm).perf_learn(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), p(learn.(tm)(sbj_i))}
            end
            searchlight.wo_anova.(version).(tm).perf_total(sbj_i,:)={sprintf('sub-%.2d',(sbj_id_list(sbj_i))), [results_wo_anova.total_perf]}
        end
        
        
    end
    
    searchlight.(version).(tm).early=cell2mat(searchlight.(version).(tm).perf_total(idx_early,2));
    searchlight.(version).(tm).late=cell2mat(searchlight.(version).(tm).perf_total(idx_late,2));
    searchlight.(version).(tm).fail=cell2mat(searchlight.(version).(tm).perf_total(idx_fail,2));
    
    searchlight.(version).(tm).mean_perf_early=mean(searchlight.(version).(tm).early);
    searchlight.(version).(tm).mean_perf_late=mean(searchlight.(version).(tm).late);
    searchlight.(version).(tm).mean_perf_fail=mean(searchlight.(version).(tm).fail);
    if t==3
        searchlight.wo_anova.(version).(tm).mean_perf_early=mean(searchlight.(version).(tm).early);
        searchlight.wo_anova.(version).(tm).mean_perf_late=mean(searchlight.(version).(tm).late);
        searchlight.wo_anova.(version).(tm).mean_perf_fail=mean(searchlight.(version).(tm).fail);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% performance plotting
    % subj별 회색 선 그리고 빨간색으로 linear fitting
    input=searchlight.(version).(tm).good_data;
    index = 1:size(input, 2);
    mean_vals = mean(input, 1);
    std_vals = std(input, 0, 1);
    
    figure;
    hold on;
    for i = 1:size(input, 1)
        plot(index, input(i, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    end
    
    fill([index, fliplr(index)], [mean_vals + std_vals, fliplr(mean_vals - std_vals)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(index, mean_vals, 'r-', 'LineWidth', 2);
    
    xlabel(tm);
    ylabel('Performance');
    title(sprintf('%s-%s-model performance',version,tm));
    
    hold off;
    %     saveas(gcf, fullfile(path_root,'OCAT_result_250101', [version, '_', tm, '.png']));
end
end

save(fullfile(path_root,'result_250102.mat'),'searchlight_trial');

%% Visualizing Brainmaps
%
% cd([path_root '\princeton-mvpa-toolbox-data-sets\tutorial_easy\working_set'])
% subj = load_afni_mask(subj,'wholebrain','wholebrain+orig');
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% plot_stability(subj,'epi','runs')
% plot_stability(subj,'epi','runs', 'proc_patname','epi_z')
%
%
%
