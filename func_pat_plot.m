% pattern similarity --> graph 그려서 저장
%% 실행 코드 예시:
% cd(~/pattern_similarity/main')
% disp(pwd); mkdir('./figures/')
% func_pat_plot(final_pat.perform)


%% function
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
                    %                     close(gcf);
                end
            end
        end
    end
end

end