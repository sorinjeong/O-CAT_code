% pattern similarity --> bilateral 계산, performance group별 나눠주고, pat structure에 저장해줌. 

%% 실행 코드 예시:
% % input 변수들
% T_L=readtable(string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet', curr_names{r});
% T_R=readtable(string(strcat((diff_same{:}),'_',(hpc_or_ctx{:}),'_',(pp{:}), '_pattern.xlsx')),'Sheet', curr_names{r+(numel(curr_names)/2)});
% pat.(pp{:}).(hpc_or_ctx{:}).(diff_same{:}).(curr_names{r})=T_L;
% pat.(pp{:}).(hpc_or_ctx{:}).(diff_same{:}).(curr_names{r+(numel(curr_names)/2)})=T_R;
% for col=1:numel(T_L(1,:))
%     T_Bi(:,col)=mean([table2array(T_L(:,col)),table2array(T_R(:,col))],2);
% end
% now.(hpc_or_ctx{:}).(pp{:}).(diff_same{:}).(bi_names{r})=T_Bi;
% T = array2table(T_Bi, 'VariableNames',T_L.Properties.VariableNames);

% % 실제 실행 코드
%[pat]=func_save_pat(T, T_L, T_R, idx_late, 'LATE', diff_same{:}, hpc_or_ctx{:}, pp{:}, bi_names, curr_names, r,pat);


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
