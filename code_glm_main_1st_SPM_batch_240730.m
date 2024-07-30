%%

% -- Modified by Sorin Jeong, July 2024

% Result Report
% Random Effect Analysis (2nd-level Model for Between-Subject Analysis
clear; clc; close all;
%set path
addpath('Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code');

%90번컴
% addpath('C:\Users\User\Documents\MATLAB\spm12')
% root_path = 'C:\Users\User\Desktop\JSR';%90번컴

%69번컴
addpath('C:\Users\leelab\Documents\MATLAB\spm12')
root_path = 'G:\JSR'; %69번컴
% addpath(genpath(root_path));%rehash path
cd(fullfile(root_path, 'spm_prep_glm_0730'))


%% Input!!

main_or_ODT = {'pre-ODT','post-ODT'}; % {'main','pre-ODT','post-ODT'} for for loops
ver='ver_1';

%% 1) Defining pathway and subjects
% load regressors
load('regressors_GLM_0730.mat');
reg=reg_0730;
sbj_id_list=sbj_id_list_38;
% set-up directory
path_in = fullfile(root_path,'spm_prep_bids', 'spmprep_final'); % fmriprep's output file (before smoothing)


%% trim_0730ming nifti and make regressor
for m_o = main_or_ODT
disp(m_o)
for i = 1:length(sbj_id_list)
    c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
    sbj_dir = fullfile(path_in, c_sbj);

    %% 1st GLM - Make multiple regressor

    mov_txt=dir(fullfile(sbj_dir,'func','rp_af*.txt'));
    total_mov=load(fullfile(mov_txt.folder,mov_txt.name));
    R_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};

    % main_or_ODT = {'main'};%,'pre-ODT','post-ODT'}

    path_nii = string(fullfile(root_path, 'spm_prep_bids','trim_0730', m_o{:}));
    if ~exist(path_nii,"dir"); mkdir(path_nii);end
    addpath(path_nii)
    current_nii_out = fullfile(path_nii,c_sbj);mkdir(current_nii_out);

    %% 4) 1st GLM - Extract 'Main' Scan
    % read 4D-scan

    % main task phase
    if strcmp(m_o{:},'main')
        scan_start = reg{1,i}.main.(ver).scan_num(1);
        scan_end = reg{1,i}.main.(ver).scan_num(2);

        valid_idx=~cellfun(@isempty, reg{1,i}.main.(ver).onset);

        names = reg{1,i}.main.(ver).name(valid_idx);
        onsets =reg{1,i}.main.(ver).onset(valid_idx);
        durations = reg{1,i}.main.(ver).duration(valid_idx);

        % pre-ODT phase
    elseif strcmp(m_o{:},'pre-ODT')
        scan_start = reg{1,i}.ODT.pre.scan_num(1);
        scan_end = reg{1,i}.ODT.pre.scan_num(2);

        valid_idx=~cellfun(@isempty, reg{1,i}.ODT.pre.regress_onset);

        names = reg{1,i}.ODT.pre.regress_name(valid_idx);
        onsets =reg{1,i}.ODT.pre.regress_onset(valid_idx);
        durations = cellfun(@(x) zeros(size(x)), onsets, 'UniformOutput', false);

        % post-ODT phase
    elseif strcmp(m_o{:},'post-ODT')
        scan_start = reg{1,i}.ODT.post.scan_num(1);
        scan_end = reg{1,i}.ODT.post.scan_num(2);
        valid_idx=~cellfun(@isempty, reg{1,i}.ODT.post.regress_onset);

        names = reg{1,i}.ODT.post.regress_name(valid_idx);
        onsets =reg{1,i}.ODT.post.regress_onset(valid_idx);
        durations = cellfun(@(x) zeros(size(x)), onsets, 'UniformOutput', false);
    end



    scan_range = scan_start:scan_end;

    files_in=arrayfun(@(x) fullfile(sbj_dir, 'func', strcat('raf', '*-*', sprintf('%.5d-%.6d',x,x), '-01.nii')), scan_range, 'UniformOutput', false);

    % create nifti file
    for s=1:length(scan_range)
        nii_data=dir(files_in{s});
        if ~isempty(nii_data)
            % nii_info.ImageSize = size(nii_data(:, :, :, scan_range));
            current_nifti = fullfile(current_nii_out,strcat(c_sbj, '_', (m_o{:}), '_',string(s), '.nii'));
            nii_info = niftiinfo(fullfile(nii_data.folder,nii_data.name));
            niftiwrite(niftiread(fullfile(nii_data.folder,nii_data.name)), current_nifti, nii_info);
        else
            scan_range=scan_range(1:end-1);

        end
    end
    if strcmp(m_o{:},'main')
        save(fullfile(path_nii, [ c_sbj, '_' (m_o{:}) '_' ver '_240730.mat']), 'names', 'onsets', 'durations');
    else
        save(fullfile(path_nii, [ c_sbj, '_' (m_o{:}) '_240730.mat']), 'names', 'onsets', 'durations');
    end
    R=total_mov(scan_range,:);

    save(fullfile(sbj_dir, ['../' c_sbj,'_' (m_o{:}) '_mov.mat']), 'R', 'R_names');


end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;
for m_o = main_or_ODT
% m_o = {'main'};
    disp(m_o)
for i = 1:length(sbj_id_list)

    c_sbj = sprintf('sub_%d', sbj_id_list(i)); disp(c_sbj)
    path_nii = fullfile(root_path, 'spm_prep_bids','trim_0730', m_o{:});
    if strcmp(m_o{:},'main')
        path_glm_out = fullfile(root_path, 'spm_prep_glm_0730', m_o{:},ver,c_sbj);
    else
        path_glm_out = fullfile(root_path, 'spm_prep_glm_0730', m_o{:},c_sbj);
    end

    if ~exist(path_glm_out,"dir"); mkdir(path_glm_out);end
    sbj_dir = fullfile(path_in, c_sbj);

    niiscans=dir(fullfile(path_nii,c_sbj, '*.nii'));
    current_nifti=cellfun(@(x) fullfile(niiscans(1).folder, x),  {niiscans.name}',  'UniformOutput',  false);

    curr_mov=fullfile(sbj_dir, ['../' c_sbj,'_' (m_o{:}) '_mov.mat']);

    if strcmp(m_o{:},'main')
        curr_reg=fullfile(path_nii, [ c_sbj, '_' (m_o{:}) '_' ver '_240730.mat']);
    else
        curr_reg=fullfile(path_nii, [ c_sbj, '_' (m_o{:}) '_240730.mat']);
    end

    %%  1st GLM - fMRI model specification
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {path_glm_out};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;

    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = current_nifti;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {curr_reg};
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


