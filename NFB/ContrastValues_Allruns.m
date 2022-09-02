%%%Script to extract contrast value (upregulation > baseline) from ROI%%%


clear all
close all
subjects= [6;7;8;9;10;11;13;14;15;17;18;20;21;22;23;24;25;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;55;57;58]; 
globalDir=pwd;
group = ['rAIC'; 'lV1'; 'rV1'; 'rAIC'; 'rAIC'; 'rAIC'; 'lV1'; 'rAIC'; 'rV1'; 'rV1'; 'lV1';'rAIC'; 'rV1';
      "rV1"; "rV1" ;"lV1";"rAIC"; "lV1";"rV1";"rAIC";"lV1";"lV1";"rAIC";"rV1";"rAIC";"lV1";"rAIC";"rV1";"rAIC";"lV1";"rV1" ; "lV1" ;"rAIC"; 
      "rV1";"rAIC" ;"rAIC" ;"rAIC"; "rAIC"; "lV1" ;"rV1";"rAIC";"rAIC";"lV1" ;"rAIC"; "rAIC"; "rAIC"];
sessDir =  {'ses-02', 'ses-03'};
statsDir = 'statsNewOn';
maskDir = 'mask'; %mask directory is inside of Session directory
con_Images = {'con_0001.nii' 'con_0002.nii' 'con_0003.nii' 'con_0004.nii' 'con_0005.nii' 'con_0006.ni' 'con_0007'};%SPM contrast image for each run
roiname = 'rvoi_rAI_mask.nii';%name of the mask, has to be in same dimension like functional files, if not normalize

for subject = 1:numel(subjects)
	subjDir = [sprintf('%03d', subjects(subject))];
    
    for session = 1:numel(sessDir)
        statsDir2 = [globalDir filesep 'sub-VP' subjDir filesep sessDir{session} filesep statsDir filesep];
        roiFilename = [globalDir filesep 'mask' filesep roiname]
        rAICmap = spm_read_vols(spm_vol(roiFilename));
        rAICvoxels=find(rAICmap); %takes the voxel number of all the nonzeros in the maskfile
        for k = 1:numel(con_Images)
            frame_names = spm_select('List', statsDir2, con_Images{k});
            frame_names3 = strcat(statsDir2,cellstr(frame_names));
            fmriTemp =  spm_read_vols(spm_vol(frame_names3{1}));
            conAc = fmriTemp(rAICvoxels);
            meanCon(k,1) = nanmean(conAc)
        end
        SessionDiff(session,:)= meanCon;
    end
    AllsubjDiff(subject,:) = [SessionDiff(1,:), SessionDiff(2,:)];
    id{subject,1} =  {['VP' subjDir]};
   
end

T = array2table(AllsubjDiff);
T.Properties.VariableNames(1:14) = {'r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11', 'r12', 'r13'};
T2 = addvars(T, id, group)
cd final_Tabel
writetable(T2,'ConnrAICSessionWholeCond.xlsx')