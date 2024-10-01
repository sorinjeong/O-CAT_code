%% activity and learning rate DTW
close all
% load('regressors_GLM_0922.mat','roi_hpc_name','roi_ctx_name')
cd('G:\JSR\240922_new_fmri\Activity')
load('G:\JSR\240922_new_fmri\bhv\data_learning_curve\curve_95/pdata_table.mat')
load('../learn_curv_info.mat')
load('activ.mat')
acq_95_41(exclud_sbj)=[];

dtw_output = struct;
dtw_dist=struct;
max_len = 0;

hpc_or_ctx={'ctx'};     %{'hpc','ctx'}

curr_names=fieldnames(activ.(hpc_or_ctx{:}));
curr_names=curr_names(~strcmp(curr_names,'trials'));

aligned_info = table;

%%
y_mat=struct;
for perform_group={'all','early','late','fail'}
    for r=1:numel(curr_names)
        
        if strcmp(perform_group,'all')
            idx=1:numel(sbj_id_list_41);
        elseif strcmp(perform_group,'early')
            idx=idx_early;
        elseif strcmp(perform_group,'late')
            idx=idx_late;
        elseif strcmp(perform_group,'fail')
            idx=idx_fail;
        end
        sbj_group=sbj_id_list_41(idx);
        sbj_learn_trial=acq_95_41(idx);
        
        
        
        % First pass to determine the maximum length
        temp_dist=[];
        for sbj_i = 1:numel(sbj_group)
            y = pdata_table.(sprintf('sub_%.2d', sbj_group(sbj_i)))';
            x = activ.(hpc_or_ctx{:}).(curr_names{r})(idx(sbj_i), :);
            
            [dist, ix, iy] = dtw(x, y);
            temp_dist = [temp_dist; dist]; % distance 값 모으기
            
            max_len = max(max_len, length(ix));
            max_len = max(max_len, length(iy));
        end
        dtw_dist.(perform_group{:}).(curr_names{r})=temp_dist;
        
        aligned_x_all = nan(numel(sbj_group), max_len);
        aligned_y_all = nan(numel(sbj_group), max_len);
        aligned_learn_trial = nan(numel(sbj_group), 1);
        
        
        for sbj_i = 1:numel(sbj_group)
            y = pdata_table.(sprintf('sub_%.2d',sbj_group(sbj_i)))';
            y_mat.(perform_group{:}).(curr_names{r})(sbj_i,:)=y;
            x=activ.(hpc_or_ctx{:}).(curr_names{r})(idx(sbj_i),:);
            
            [dist, ix, iy] = dtw(x, y);
            
            dtw_output.(perform_group{:}).(curr_names{r}).(sprintf('sub_%.2d', sbj_group(sbj_i))) = {'dist', dist; 'ix', ix'; 'iy', iy'};
            
            aligned_x_all(sbj_i, 1:length(ix)) = x(ix);
            aligned_y_all(sbj_i, 1:length(iy)) = y(iy);
            
            if ~isnan(sbj_learn_trial(sbj_i))
                aligned_learn_trial(sbj_i) = find(iy==sbj_learn_trial(sbj_i),1);
            end
        end
        
        mean_aligned_x.(perform_group{:}).(curr_names{r}) = nanmean(aligned_x_all, 1);
        mean_aligned_y.(perform_group{:}).(curr_names{r}) = nanmean(aligned_y_all, 1);
        if ~isempty(aligned_learn_trial)
            mean_trial.(perform_group{:}).(curr_names{r}) = nanmean(aligned_learn_trial, 1);
        end
    end
end

%% plotting
%% plotting original signals

plot_count = 0;
figure_index = 1;
for perform_group={'all','early','late','fail'}
    for r=1:numel(curr_names)
        
        
        [h_e_l,p_e_l]=ttest2(dtw_dist.early.(curr_names{r}), dtw_dist.late.(curr_names{r}));
        fprintf('Early vs Late DTW distance t-test: h = %d, p = %.4f\n', h_e_l, p_e_l);
        [h_e_f,p_e_f]=ttest2(dtw_dist.early.(curr_names{r}), dtw_dist.fail.(curr_names{r}));
        fprintf('Early vs Fail DTW distance t-test: h = %d, p = %.4f\n', h_e_f, p_e_f);
        [h_l_f,p_l_f]=ttest2(dtw_dist.late.(curr_names{r}), dtw_dist.fail.(curr_names{r}));
        fprintf('Late vs Fail DTW distance t-test: h = %d, p = %.4f\n', h_l_f, p_l_f);
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            fprintf('p-valueeeeeeeeeeeeee!!!!!!!!!!!!!!!!\n')
        end
        
        if mod(plot_count, 6) == 0
            figure('position',[827,536,1390,725]);
            figure_index = figure_index + 1;
        end
        
        plot_count = plot_count + 1;
        
        % Plot original signals
        subplot(3, 2, mod(plot_count-1, 6) + 1);
        yyaxis left;
        plot(mean(activ.(hpc_or_ctx{:}).(curr_names{r}), 1), '-o');
        ylabel('activity');
        
        yyaxis right;
        plot(mean(y_mat.(perform_group{:}).(curr_names{r}), 1), '-x');
        ylabel('learning rate');
        title(sprintf('Mean Original Signals (%s - %s)', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.')));
        
        if mod(plot_count, 3) == 0
            legend('activity', 'learning rate','Position',[0.889172287757767,0.015557516378381,0.094834884150769,0.050065018188194]);
        elseif mod(plot_count, 6) == 4
            legend('activity', 'learning rate','Position',[0.887014014376472,0.611419585343898,0.094834884150769,0.050065018188194]);
        end
        
        %% plotting aligned signals
        if mod(plot_count, 6) == 0
            figure('position',[827,536,1390,725]);
            figure_index = figure_index + 1;
        end
        
        plot_count = plot_count + 1;
        
        subplot(3, 2, mod(plot_count-1, 6) + 1);
        yyaxis left;
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            plot(mean_aligned_x.(perform_group{:}).(curr_names{r}), '-o', 'LineWidth', 2); % p value가 0.1 이하일 때 선 굵게 표시
        else
            plot(mean_aligned_x.(perform_group{:}).(curr_names{r}), '-o');
        end
        ylabel('activity (aligned)');
        
        yyaxis right;
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            plot(mean_aligned_y.(perform_group{:}).(curr_names{r}), '-x', 'LineWidth', 2); % p value가 0.1 이하일 때 선 굵게 표시
        else
            plot(mean_aligned_y.(perform_group{:}).(curr_names{r}), '-x');
        end
        ylabel('learning rate (aligned)');
        % Add additional information to the title
        additional_info = '';
        if p_e_l < 0.1
            additional_info = [additional_info, 'Early vs Late, '];
        end
        if p_e_f < 0.1
            additional_info = [additional_info, 'Early vs Fail, '];
        end
        if p_l_f < 0.1
            additional_info = [additional_info, 'Late vs Fail, '];
        end
        if ~isempty(additional_info)
            additional_info = [' (', additional_info(1:end-2), ')'];
        end
        
        title(sprintf('Mean Aligned Signals (%s - %s)%s', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.'), additional_info));
        
        % mean_learning trial을 세로로 그리기
        if exist('mean_trial', 'var') && ~isnan(mean_trial.(perform_group{:}).(curr_names{r})) && ~isempty(mean_trial.(perform_group{:}).(curr_names{r}))
            xline(mean_trial.(perform_group{:}).(curr_names{r}), '--k', 'LineWidth', 1.5);
        end
        
        
        save_folder=sprintf('G:/JSR/240922_new_fmri/Activity/DTW/%s_sub',perform_group{:});
        mkdir(save_folder)
        saveas(gcf, sprintf('%s/DTW_%s_%d.png', save_folder,hpc_or_ctx{:},figure_index-1));
        
    end
end


%% Smoothen and normalized neural activity

% 함수 정의
%% 정규화 및 스무딩 함수
normalize = @(x) (x - min(x)) / (max(x) - min(x));
smooth = @(x) movmean(x, 5); % 이동 평균을 사용한 스무딩 (5씩 이동)


plot_count = 0;
figure_index = 1;
for perform_group={'all','early','late','fail'}
    for r=1:numel(curr_names)
        
        [h_e_l,p_e_l]=ttest2(dtw_dist.early.(curr_names{r}), dtw_dist.late.(curr_names{r}));
        fprintf('Early vs Late DTW distance t-test: h = %d, p = %.4f\n', h_e_l, p_e_l);
        [h_e_f,p_e_f]=ttest2(dtw_dist.early.(curr_names{r}), dtw_dist.fail.(curr_names{r}));
        fprintf('Early vs Fail DTW distance t-test: h = %d, p = %.4f\n', h_e_f, p_e_f);
        [h_l_f,p_l_f]=ttest2(dtw_dist.late.(curr_names{r}), dtw_dist.fail.(curr_names{r}));
        fprintf('Late vs Fail DTW distance t-test: h = %d, p = %.4f\n', h_l_f, p_l_f);
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            fprintf('p-valueeeeeeeeeeeeee!!!!!!!!!!!!!!!!\n')
        end
        
        if mod(plot_count, 6) == 0
            figure('position',[827,536,1390,725]);
            figure_index = figure_index + 1;
        end
        
        plot_count = plot_count + 1;
        
        % Plot original signals
        subplot(3, 2, mod(plot_count-1, 6) + 1);
        yyaxis left;
        plot(smooth(normalize(mean(activ.(hpc_or_ctx{:}).(curr_names{r}), 1))), '-o');
        ylabel('activity');
        
        yyaxis right;
        plot(smooth(normalize(mean(y_mat.(perform_group{:}).(curr_names{r}), 1))), '-x');
        ylabel('learning rate');
        title(sprintf('Mean Original Signals (%s - %s)', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.')));
        
        if mod(plot_count, 3) == 0
            legend('activity', 'learning rate','Position',[0.889172287757767,0.015557516378381,0.094834884150769,0.050065018188194]);
        elseif mod(plot_count, 6) == 4
            legend('activity', 'learning rate','Position',[0.887014014376472,0.611419585343898,0.094834884150769,0.050065018188194]);
        end
        
        %% plotting aligned signals
        if mod(plot_count, 6) == 0
            figure('position',[827,536,1390,725]);
            figure_index = figure_index + 1;
        end
        
        plot_count = plot_count + 1;
        
        subplot(3, 2, mod(plot_count-1, 6) + 1);
        yyaxis left;
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            plot(smooth(normalize(mean_aligned_x.(perform_group{:}).(curr_names{r}))), '-o', 'LineWidth', 2);
        else
            plot(smooth(normalize(mean_aligned_x.(perform_group{:}).(curr_names{r}))), '-o');
        end
        ylabel('activity (aligned)');
        
        yyaxis right;
        if p_e_l < 0.1 || p_e_f < 0.1 || p_l_f < 0.1
            plot(smooth(normalize(mean_aligned_y.(perform_group{:}).(curr_names{r}))), '-x', 'LineWidth', 2);
        else
            plot(smooth(normalize(mean_aligned_y.(perform_group{:}).(curr_names{r}))), '-x');
        end
        ylabel('learning rate (aligned)');
        
        % Add additional information to the title
        additional_info = '';
        if p_e_l < 0.1
            additional_info = [additional_info, 'Early vs Late, '];
        end
        if p_e_f < 0.1
            additional_info = [additional_info, 'Early vs Fail, '];
        end
        if p_l_f < 0.1
            additional_info = [additional_info, 'Late vs Fail, '];
        end
        if ~isempty(additional_info)
            additional_info = [' (', additional_info(1:end-2), ')'];
        end
        
        title(sprintf('Mean Aligned Signals (%s - %s)%s', hpc_or_ctx{:}, strrep(curr_names{r}, '_', '.'), additional_info));
        % mean_learning trial을 세로로 그리기
        if exist('mean_trial', 'var') && ~isnan(mean_trial.(perform_group{:}).(curr_names{r})) && ~isempty(mean_trial.(perform_group{:}).(curr_names{r}))
            xline(mean_trial.(perform_group{:}).(curr_names{r}), '--k', 'LineWidth', 1.5);
        end
        
        
        save_folder=sprintf('G:/JSR/240922_new_fmri/Activity/DTW/smoothing/%s_sub',perform_group{:});
        mkdir(save_folder)
        saveas(gcf, sprintf('%s/DTW_%s_%d.png', save_folder,hpc_or_ctx{:},figure_index-1));
        
    end
end


% aligned_info.Properties.VariableNames = {'hpc_cortex', 'ROI', 'perform_group', 'trial_neural_activity', 'trial_learn'};
%
% aligned_info
%
% writetable(aligned_info,'./DTW/aligned_info.xlsx')