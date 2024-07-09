# O-CAT_code

# OCAT - Behavior
(last modified : 2024.07.09)

**1. code_bhv_logparsing_func.m**
  - func_bhv_logparsing.m
  - func_create_temp_table.m
  - func_perform_graph.m
  - func_bhv_performance_plots.m

output:
- data_bhv_log_table
  - (sub-**)_task-ocat_run-01_events.tsv
  - event_table
  - event_TR
  - event_pre_PV
  - event_post_PV
  - num_sbj_events

- data_bhv_plot
  - accuracy plot
  - RT plot
  - bias plot
  - performance plot
   
 - data_learning_curve
   - responses as binary array for learning curve input  

**2. regressors_GLM_0709.m**
   : make glm regressors from (sub-**)_task-ocat_run-01_events.tsv files

output:
- all_subjects_task-ocat_run-01_events.xlsx
  - added [NavEnd, ITIEnd]
- regressors_GLM_0709.mat
  - "reg_for_glm_ver1","reg_for_glm_ver2","sbj_id_list_38"
 
**3. code_glm_main_1st_SPM_batch_240619.m**





