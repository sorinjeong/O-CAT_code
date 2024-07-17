% batch_newLSS

% India: Single trial regressor
% April, 2023

% General settings
% tStart = tic;
% load subject list
% load('D:\Delta_T_Analysis\Behavioral\Subject_data_final\Final_sbj_run_list.mat');
% sbj_run_list(2,:) = [];
% subjects = sbj_run_list;
% subject = {'709'};
% subject = '709';

% spmFolder = '/nfs/ep2/AX/first_levels/00_MONTH/';
% spmSubFolder = '/MTU/func_an_SPM8/';
%
% outFolder = '/home/tsalo/lssForUnivariate/';
% outSubFolder = '/';

% LSS settings
% includeConditions = {'Asso' 'TOJ_close' 'TOJ_far'};
% settings.model = 2;            % 1- Rissman, 2- LSS
% settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
% settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% LSS settings
% includeConditions = {'CueA' 'CueB'};
% settings.model = 2;            % 1- Rissman, 2- LSS
% settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
% settings.overwrite = 0;        % 0- do not overwrite, 1- overwrite

% Connectivity settings (not used in my analysis)
% rois = load('/nfs/cntracs/lssForKim/gt_rois_final.mat');
% settings.fConnType = 'roi2roi'; % seed2voxel or roi2roi
%
% for iSubj = 1:length(subjects)
%     spmDir = [spmFolder subjects{iSubj} spmSubFolder];
%     outDir = [outFolder subjects{iSubj} outSubFolder];
%
%     images = lssGenerateBetasSpm(subjects{iSubj}, spmDir, outDir, includeConditions, settings);
%     lssCorrelation(images, rois.rois, settings);
% end


%% Set path & execute lssGenerateBetasSPM function
% for extracting parameter estimates from each ROIs for each trials
% for extracting patterns from each ROIs for each trials


%% Run lssGenerateBetasSPM.m
clear; clc; close all;
%% path
% delete(gcp('nocreate'))
% restoredefaultpath;
% rehash toolboxcache;
root_path = 'C:\Users\User\Desktop\JSR';%90번컴
addpath('C:\Users\User\Documents\MATLAB\spm12')
addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code')
cd(root_path)
%% Get data

% regressor=reg_0711
load('regressors_GLM_0717.mat','sbj_id_list_38');
sbj_id_list=sbj_id_list_38;


   settings.model = 2;            % 1- Rissman, 2- LSS
           settings.useTempFS = 0;        % 0- do not use temporary files, 1- use temporary files
                   settings.overwrite = 1;        % 0- do not overwrite, 1- overwrite


%%
% parfor main_or_ODT = {'main','pre-ODT','post-ODT'}
    % main_or_ODT = {'post-ODT'};
main_or_ODT = {'main'}

    path_out =  fullfile(root_path, 'spm_prep_bids','trim',main_or_ODT{:});

    regressor=load(fullfile(path_out, [ sprintf('sub_%d', 1), '_' (main_or_ODT{:}) '_regressor_240717.mat']));

    includeConditions = regressor.names(1:4);

    for sbj_i = 1:length(sbj_id_list)
        % sbj_i = 1
        c_sbj = sprintf('sub_%d', sbj_id_list(sbj_i)); disp(c_sbj)
        subject=c_sbj;
        %         path_out =  fullfile(root_path, 'spm_prep_glm',main_or_ODT{:});
        %
        %         load(fullfile(path_out, [ c_sbj, '_' (main_or_ODT{:}) '_regressor_240711.mat']));



        % change this for different Design Matrix
        spmDir = fullfile(root_path, 'spm_prep_glm', main_or_ODT{:},c_sbj);
        outDir = fullfile(root_path, 'spm_prep_glm', main_or_ODT{:}, 'single_trial', c_sbj);

        % spmDir = ['C:\Users\zyeon\SynologyDrive\OCAT\Results_new\GLM\1stLevel\glm' num2str(glm_i) '\SUB' num2str(sbj_id_list(sbj_i))];
        % outDir = ['C:\Users\zyeon\SynologyDrive\OCAT\Results_new\GLM\1stLevel\SinggleTrial\glm' num2str(glm_i) '\SUB' num2str(sbj_id_list(sbj_i))];

        % LSS settings
        % includeConditions = behav_regressor_glm{glm_i}{sbj_i}.regress_name;
        %         includeConditions = includeConditions(ismember(includeConditions, oc_list_name));
        %         includeConditions = names;


          % myCluster = parcluster('Processes');
        % delete(myCluster.Jobs)
        delete(gcp('nocreate'));

        lssGenerateBetasSpm(c_sbj, spmDir, outDir, includeConditions, settings)

    end  %-- end of for ct_sub





