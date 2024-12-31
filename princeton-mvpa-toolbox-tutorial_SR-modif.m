close all; clc; clear
path_root='G:\JSR\241104_new_fmri\activity\princeton_searchlight';
addpath(genpath(path_root))
cd([path_root '\princeton-mvpa-toolbox-master'])

mvpa_add_paths;
help tutorial_easy 

%'./core/preproc/'에 있음
mex './core/preproc/compute_xcorr.c'

load('tutorial_regs.mat')
load('tutorial_runs.mat')


% open './core/learn/perfmet_maxclass.m'
% open './core/template/perfmet_template.m'

%% 1) Initializing the subject (''subj'') structure
% open('init_subj.m') 한 후 datetime(true) --> datetime('now')로 변경!!
subj = init_subj('haxby8','tutorial_subj');
summarize(subj);

%% 2) Ventral temporal with GLM mask
% 1) open('init_object.m') 한 후 datetime(true) --> datetime('now')로 변경!!
% 2) open('set_mat.m') 한 후 datetime(true) --> datetime('now')로 변경!!

% AFNI버전
% cd([path_root '\princeton-mvpa-toolbox-data-sets\tutorial_easy\working_set'])
% subj = load_afni_mask(subj,'VT_category-selective','mask_cat_select_vt+orig');
% SPM버전
cd([path_root '\princeton-mvpa-toolbox-data-sets\tutorial_easy_spm\working_set'])
subj = load_spm_mask(subj,'VT_category-selective','mask_cat_select_vt.nii');

summarize(subj)
mask = get_object(subj,'mask','VT_category-selective') 


%% 3) EPI pattern

% % AFNI버전
% for i=1:10 
% raw_filenames{i} = sprintf('haxby8_r%i+orig',i);
% end
% subj = load_afni_pattern(subj,'epi','VT_category-selective',raw_filenames);
% summarize(subj)

% SPM버전
for i=1:10  
     index=num2str(i);
     raw_filenames{i} = ['haxby8_r' index '.img'];
end
subj = load_spm_pattern(subj,'epi','VT_category-selective',raw_filenames);

pats = get_object(subj,'pattern','epi') 


%% 4) Condition regressors
cd([path_root '\princeton-mvpa-toolbox-master'])

subj = init_object(subj,'regressors','conds');
load('tutorial_regs');
whos regs

subj = set_mat(subj,'regressors','conds',regs);

condnames = {'face','house','cat','bottle','scissors','shoe','chair','scramble'};
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

condsobj= get_object(subj,'regressors','conds') 


%% 5) Runs selector
subj = init_object(subj,'selector','runs'); 
load('tutorial_runs'); 

subj = set_mat(subj,'selector','runs',runs); 
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

class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;

[subj results] = cross_validation(subj,'epi_z','conds','runs_xval','epi_z_thresh0.05',class_args);

%% 결과 확인
%total_perf 확인 (iteration 전체의 평균값)
results
% acts, perfmet 확인
results.iterations(1)
%guesses에서 예측한 결과를 확인할 수 있고, 그 퍼포먼스는 아래 코드로 확인 가능
results.iterations(1).perfmet.perf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualizing Brainmaps

cd([path_root '\princeton-mvpa-toolbox-data-sets\tutorial_easy\working_set'])
subj = load_afni_mask(subj,'wholebrain','wholebrain+orig');

















plot_stability(subj,'epi','runs')
plot_stability(subj,'epi','runs', 'proc_patname','epi_z')



