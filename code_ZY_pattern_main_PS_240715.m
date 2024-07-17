%% single trial naming
zcode_path='Z:\E-Phys Analysis\fMRI_ocat\OCAT_DIR\code';
addpath(zcode_path)
load('num_sbj_events.mat');
load('regressors_GLM_0714.mat')
sbj_id_list=sbj_id_list_38;
singletrial_path='C:\Users\User\Desktop\JSR\spm_prep_glm\GLM\main\single_trial';

roi_selection=21:26;


%% rename_ trial


% for sbj_i=1:numel(sbj_id_list)
%     sbj_n=sprintf('sub_%d',sbj_id_list(sbj_i));
%     sbj_OC=dir(fullfile(singletrial_path,sbj_n,'betas','Sess001','obj_show_*'));
%     roc=cellfun(@(x) extractBetween(x,'obj_show_','_0'),{sbj_OC.name});
%     d=reg_0714{1, sbj_i}.trial_detail;
%     un=unique(roc);
%     sidx=[];
%     for r=1:numel(un)
%         idx=d.(sprintf('trial_%s',un{r}));
%         sidx=[sidx idx];
%     end
%     [b,iii]=sort(sidx);
%     if max(iii)==32
%     tbl.roc(tbl.session==sbj_id_list(sbj_i))=roc(iii);
%     for i=1:height(sbj_OC)
%         beta=fullfile(sbj_OC(i).folder,sbj_OC(i).name,'beta_0001.nii');
%         beta_rename=sprintf('%.2d_%s.nii',sidx(i), sbj_OC(i).name);
%         copyfile(beta, fullfile(singletrial_path,sbj_n, beta_rename));
%     end
%     else
%         fprintf('error')
%     end
% end






%% ODT - Trial pattern extraction
% % trial별 roc name (hit, miss, corr_rej, false)
% load("sbj_events.mat")
%
% roc_cell=cell(1, ceil(length(roc_trial)/32));
% for i = 1:length(roc_cell)
%     startIdx = (i-1)*32 + 1;
%     endIdx = min(i*32, length(roc_trial));
%     roc_cell{i} = roc_trial(startIdx:endIdx);
% end

% % trial별 object 및 context 정보
% obj_num=cell(1, ceil(height(tbl)/32));
% context=cell(1, ceil(height(tbl)/32));
% for i = 1:length(obj_num)
%     startIdx = (i-1)*32 + 1;
%     endIdx = min(i*32, height(tbl));
%     obj_num{i} = tbl.Object(startIdx:endIdx)';
%     for cont=1:32
%         if tbl.Association(cont) ==1
%             ctxt{cont} = tbl.Context(cont);
%         elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==1
%             ctxt{cont} = 2;
%         elseif tbl.Association(cont) ==0 && tbl.Context(cont) ==2
%             ctxt{cont} = 1;
%         end
%         context{i}=ctxt;
%     end
% end
% for sbj_i=1:numel(sbj_id_list)
%     reg_0714{1, sbj_i}.trial_detail.context=cell2mat(context{sbj_i});
%     reg_0714{1, sbj_i}.trial_detail.object=obj_num{sbj_i};
% end


%% pattern structure
% load('roi_name.mat')
roi_name=seg.seg_fit.hpc.roi_name_list;
curr_names=roi_name(roi_selection);

OC_pattern=struct;

for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    load(strcat('./ZY_segmentation/data_fmri_organized_seg_v5/', string(sbj_n)))

    roi_nn=seg.seg_fit.hpc.roi_list_nn;

    sbj_nii=dir(fullfile(singletrial_path,sprintf('sub_%d',sbj_n),'*.nii'));
    for r=1:numel(curr_names)
        for t=1:height(sbj_nii) % 32 trials
            curr_beta=niftiread(fullfile(sbj_nii(1).folder,sbj_nii(t).name));
            curr_mask=roi_nn{roi_selection(r)};
            pattern_roi=curr_beta(curr_mask);

            OC_pattern.data{1,sbj_i}{1,r}{1,t}=pattern_roi;
            OC_pattern.data{1,sbj_i}{1,r}{2,t}=extractBetween(sbj_nii(t).name,'obj_show_','_0');
            OC_pattern.data{1,sbj_i}{1,r}{3,t}=obj_num{1,sbj_n}(t);%object정보
            OC_pattern.data{1,sbj_i}{1,r}{4,t}=context{1,sbj_n}{t}; %context정보
            OC_pattern.data{1,sbj_i}{2,r}=roi_name{roi_selection(r)};
        end

    end
    OC_pattern.raw_roi=roi_nn;
end
% subject(27) - roi region(30) - trial(32)


save('C:\Users\User\Desktop\JSR\spm_prep_glm\OC_pattern.mat',"OC_pattern");


%%%%%%%%%%%%%%%%%%%% 이 아래 OC pattern에 맞춰서 다시 수정하기!!!
%% same Ps, diff PS
parsing=struct;load('C:\Users\User\Desktop\JSR\spm_prep_glm\OC_pattern.mat');
% load('PS_basic_info.mat');
sbj_id_list=sbj_id_list(1:28)
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    for r=1:numel(roi_selection)

        curr_s_r=OC_pattern.data{1,sbj_i}{1,r};
        bulk_corr=corrcoef(cell2mat(curr_s_r(1,:)), 'rows','pairwise');


        f_idx = find(cellfun(@(x) x == 1, curr_s_r(4,1:4)));
        curr_f= cell2mat(OC_pattern.data{1,sbj_i}{1,r}(3,f_idx));
        curr_c= setdiff(4:7,curr_f);

        f1_idx=find(cell2mat(curr_s_r(3,:))==curr_f(1));
        f2_idx=find(cell2mat(curr_s_r(3,:))==curr_f(2));
        c1_idx=find(cell2mat(curr_s_r(3,:))==curr_c(1));
        c2_idx=find(cell2mat(curr_s_r(3,:))==curr_c(2));

        %% same
        same_ctxt.all{r,sbj_i}=[bulk_corr(f1_idx,f2_idx); bulk_corr(c1_idx,c2_idx)];

        m=[];
        for lap=1:8
            m(lap)=mean([same_ctxt.all{r,sbj_i}(lap,lap); same_ctxt.all{r,sbj_i}(lap+8,lap)],"all");
        end
        same_ctxt.lap{1,sbj_i}{:,r}=m';
        for block=1:4
            same_ctxt.block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
        end

        %% diff
        % complex 1: f1 - c1 / complex 2: f1 - c2 / complex 3: f2-c1 / complex 4: f2-c2

        diff_ctxt.all{r,sbj_i}{1}=bulk_corr(f1_idx,c1_idx);
        diff_ctxt.all{r,sbj_i}{2}=bulk_corr(f1_idx,c2_idx);
        diff_ctxt.all{r,sbj_i}{3}=bulk_corr(f2_idx,c1_idx);
        diff_ctxt.all{r,sbj_i}{4}=bulk_corr(f2_idx,c2_idx);

        m=[];

        for lap=1:8
            m(lap)=mean([diff_ctxt.all{r,sbj_i}{1}(lap,lap); diff_ctxt.all{r,sbj_i}{2}(lap,lap);diff_ctxt.all{r,sbj_i}{3}(lap,lap);diff_ctxt.all{r,sbj_i}{4}(lap,lap)],"all");
        end
        diff_ctxt.lap{1,sbj_i}{:,r}=m';
        for block=1:4
            diff_ctxt.block{1,sbj_i}{block,r}=mean(m(block*2-1:block*2));
        end

    end
    same_ctxt.half{1,sbj_i}(1,:)=mean(cell2mat(same_ctxt.block{1,sbj_i}(1:2,:)));
    same_ctxt.half{1,sbj_i}(2,:)=mean(cell2mat(same_ctxt.block{1,sbj_i}(3:4,:)));
    diff_ctxt.half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(1:2,:)));
    diff_ctxt.half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(3:4,:)));
end

save('patterns_for_analysis_0715',"patterns","same_ctxt","diff_ctxt","sbj_id_list")

roi_name(roi_selection)
%%
m2=[];            m2=mean([mean(diff_ctxt.all{r,sbj_i}{1}); mean(diff_ctxt.all{r,sbj_i}{2}); mean(diff_ctxt.all{r,sbj_i}{3}); mean(diff_ctxt.all{r,sbj_i}{4})]);
new_mean.lap{1,sbj_i}{:,r}=m2;
 for block=1:4
            new_mean.block{1,sbj_i}{block,r}=mean(m2(block*2-1:block*2));
 end
 new_mean.half{1,sbj_i}(1,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(1:2,:)));
 new_mean.half{1,sbj_i}(2,:)=mean(cell2mat(diff_ctxt.block{1,sbj_i}(3:4,:)));


%%

now=struct;
for sbj_i=1:numel(sbj_id_list)
    sbj_n = sbj_id_list(sbj_i);
    for r=1:numel(roi_selection)

        now.diff.(roi_name{roi_selection(r)}){sbj_i,:}=diff_ctxt.half{1,sbj_i}(:,r)';
        now.same.(roi_name{roi_selection(r)}){sbj_i,:}=same_ctxt.half{1,sbj_i}(:,r)';
    end
end





%
% % --> 이거로 69번 컴에서 graph pad prism 으로 plot그리기 (block별, half 별)
%%