%%

% -- Modified by Sorin Jeong, July 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');
addpath('C:\Users\User\Documents\MATLAB\spm12')
root_path = 'C:\Users\User\Desktop\JSR';%90번컴
addpath(genpath(root_path));rehash path
cd(root_path)
%% Input!!

% main_or_ODT = 'pre-ODT'; % if you want to analize MAIN task, put 'main' or for ODT, put 'pre-ODT' or 'post-ODT'


%% 1) Defining pathway and subjects
% load regressors
load('regressors_GLM_0711.mat');
reg=reg_0711;
sbj_id_list=sbj_id_list_38;
% set-up directory
path_in = fullfile(root_path,'spm_prep_bids', 'derivatives'); % fmriprep's output file (before smoothing)


%% trimming nifti and make regressor
for i = 1:length(sbj_id_list)
    c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
    sbj_dir = fullfile(path_in, c_sbj);

    mov_txt=dir(fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','rp_af*.txt'));
    total_mov=load(fullfile(mov_txt.folder,mov_txt.name));
    R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};

    for main_or_ODT = {'main','pre-ODT','post-ODT'}
        path_out = string(fullfile(root_path, 'spm_prep_glm', main_or_ODT{:}));
        if ~exist(path_out,"dir"); mkdir(path_out);end
        addpath(path_out)
        current_beta_out = fullfile(path_out,c_sbj);mkdir(current_beta_out);

        %% 4) 1st GLM - Extract 'Main' Scan
        % read 4D-scan

        % main task phase
        if strcmp(main_or_ODT{:},'main')
            scan_start = reg{1,i}.main.scan_num(1);
            scan_end = reg{1,i}.main.scan_num(2);

            valid_idx=~cellfun(@isempty, reg{1,i}.main.regress_onset);

            names = reg{1,i}.main.regress_name(valid_idx);
            onsets =reg{1,i}.main.regress_onset(valid_idx);
            durations = cellfun(@(x) zeros(size(x)), onsets, 'UniformOutput', false);

            % pre-ODT phase
        elseif strcmp(main_or_ODT{:},'pre-ODT')
            scan_start = reg{1,i}.ODT.pre.scan_num(1);
            scan_end = reg{1,i}.ODT.pre.scan_num(2);

            valid_idx=~cellfun(@isempty, reg{1,i}.ODT.pre.regress_onset);

            names = reg{1,i}.ODT.pre.regress_name(valid_idx);
            onsets =reg{1,i}.ODT.pre.regress_onset(valid_idx);
            durations = cellfun(@(x) zeros(size(x)), onsets, 'UniformOutput', false);

            % post-ODT phase
        elseif strcmp(main_or_ODT{:},'post-ODT')
            scan_start = reg{1,i}.ODT.post.scan_num(1);
            scan_end = reg{1,i}.ODT.post.scan_num(2);

            valid_idx=~cellfun(@isempty, reg{1,i}.ODT.post.regress_onset);

            names = reg{1,i}.ODT.post.regress_name(valid_idx);
            onsets =reg{1,i}.ODT.post.regress_onset(valid_idx);
            durations = cellfun(@(x) zeros(size(x)), onsets, 'UniformOutput', false);
        end

        scan_range = scan_start:scan_end;

        files_in=arrayfun(@(x) fullfile(root_path, 'spm_prep_bids','spmprep_v5', c_sbj, 'func', strcat('raf', '*-*', sprintf('%.5d-%.6d',x,x), '-01.nii')), scan_range, 'UniformOutput', false);

        % create nifti file
        for s=1:length(scan_range)
            nii_data=dir(files_in{s});
            if ~isempty(nii_data)
                % nii_info.ImageSize = size(nii_data(:, :, :, scan_range));
                current_nifti = fullfile(current_beta_out,strcat(c_sbj, '_', (main_or_ODT{:}), '_',string(s), '.nii'));
                nii_info = niftiinfo(fullfile(nii_data.folder,nii_data.name));
                niftiwrite(niftiread(fullfile(nii_data.folder,nii_data.name)), current_nifti, nii_info);
            else
                scan_range=scan_range(1:end-1);

            end
        end
        save(fullfile(path_out, [ c_sbj, '_' (main_or_ODT{:}) '_regressor_240711.mat']), 'names', 'onsets', 'durations');

        R=total_mov(scan_range,:);

        save(fullfile(sbj_dir, [c_sbj,'_' (main_or_ODT{:}) '_mov.mat']), 'R', 'R_names');


    end
end

%
%
%% 1st GLM - Make multiple regressor
% for i = 1:length(sbj_id_list)
%     c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
%     sbj_dir = fullfile(path_in, c_sbj);
%
%     mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','rp_af*.txt');
%     total_mov=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%
%     save(fullfile(sbj_dir, [c_sbj, '_movements.mat']), 'R', 'R_names');
%
% mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','pre_mov.txt');
%     R=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_pre_mov.mat']), 'R', 'R_names');
%
% mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','post_mov.txt');
%     R=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_post_mov.mat']), 'R', 'R_names');
% end

% for i = 1:length(sbj_id_list)
%     c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
%     sbj_dir = fullfile(path_in, c_sbj);
%
%     mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func',sprintf('%d_movement_file.txt',sbj_id_list(i)));
%     R=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_movements.mat']), 'R', 'R_names');
%
% mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','pre_mov.txt');
%     R=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_pre_mov.mat']), 'R', 'R_names');
%
% mov_txt=fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','post_mov.txt');
%     R=load(mov_txt);
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_post_mov.mat']), 'R', 'R_names');
% end

% for i = 1:length(sbj_id_list)
%     c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
%     sbj_dir = fullfile(path_in, c_sbj);
% delete(fullfile(sbj_dir, [c_sbj, '_movements.mat']))
% delete(fullfile(sbj_dir, [c_sbj, '_pre_mov.mat']))
% delete(fullfile(sbj_dir, [c_sbj, '_post_mov.mat']))
% end

% for i = 29:length(sbj_id_list)
%     c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
%     sbj_dir = fullfile(path_in, c_sbj);
%
%     old=dir(fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','*mov*.txt'));
%     cellfun(@(x) delete(fullfile(old(1).folder, x)),{old.name})
%
%     mov_txt=dir(fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func','rp*-01.txt'));
%
%     R=load(fullfile(mov_txt.folder,mov_txt.name));
%     writematrix(R,fullfile(root_path,'spm_prep_bids','spmprep_v5' ,c_sbj,'func',sprintf('%d_movement_file.txt',sbj_id_list(i))),"Delimiter",';');
%
%     R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
%     save(fullfile(sbj_dir, [c_sbj, '_movements.mat']), 'R', 'R_names');
%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

for i = 1:length(sbj_id_list)
    for main_or_ODT = {'main','pre-ODT','post-ODT'}

        c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
        path_out=fullfile(root_path, 'spm_prep_glm', main_or_ODT{:});
        path_glm_out = fullfile(root_path, 'spm_prep_glm', 'GLM',main_or_ODT{:},c_sbj);
        if ~exist(path_glm_out,"dir"); mkdir(path_glm_out);end
        sbj_dir = fullfile(path_in, c_sbj);

% temp=dir(fullfile(path_out,c_sbj, '*trimmed.nii'));
% cellfun(@(x) delete(fullfile(temp(1).folder,x)),{temp.name});

        niiscans=dir(fullfile(path_out,c_sbj, '*.nii'));
        current_nifti=cellfun(@(x) fullfile(niiscans(1).folder, x),  {niiscans.name}',  'UniformOutput',  false);

        curr_mov=fullfile(sbj_dir, [c_sbj,'_' (main_or_ODT{:}) '_mov.mat']);

        %%  1st GLM - fMRI model specification
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_spec.dir = {path_glm_out};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;

        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = current_nifti;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(path_out, [ c_sbj, '_' (main_or_ODT{:}) '_regressor_240711.mat'])};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {curr_mov};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

        % load(fullfile(sbj_dir, [c_sbj, '_movements.mat']))

        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1; %Model Interactions (Volterra) / 1==no interaction , 2==interaction
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; %autoregrssive ; AR(1) is the default option (by AHN)
        spm_jobman('run', matlabbatch)


        %% 7) 1st GLM - Model estimation
        matlabbatch = [];
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(path_glm_out, 'SPM.mat') };
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch)


    end

end


