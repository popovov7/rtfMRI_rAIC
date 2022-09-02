%%This script calculates from beta images the difference between
%%Upregulation and Basline for the Baslinerun from Sess-02 and Transferrun
%%from Sess-03
%%updated 01-08-2022

clear all
subjects= [6;7;8;9;10;11;12;13;14;15;16;17;18;20;21;22;23;24;25;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;55;57;58];
globalDir = pwd;
sessDir = {'ses-02' 'ses-03'};%First and second NFB session
statsDir = 'statsNewOn';%SPM First_Level Stats Output
maskDir = 'mask'; %mask directory 
beta_images_B = {'beta_0001.nii' 'beta_0002'};%Regressors corresponding to Baseline run
beta_images_T = { 'beta_0055.nii' 'beta_0056'};%Regressors corresponding to Transfer run
roiname = 'rvoi_rAI_mask.nii';%'rvoi_V1_mask.nii


for subject = 1:numel(subjects)
        subjDir = [sprintf('%03d', subjects(subject))];
        root_f = [globalDir filesep 'sub-VP' subjDir filesep sessDir{1} filesep statsDir filesep];%First NFB session
        roiFilename = [globalDir filesep 'mask' filesep roiname];
        rAICmap = spm_read_vols(spm_vol(roiFilename));
        rAICvoxels=find(rAICmap); %takes the voxel number of all the nonzeros in the maskfile
        
        %%Baseline run
        %Baseline regressor
        frame_names = spm_select('List', root_f, beta_images_B{1});%selects the baseline regressor for baseline run
        frame_names3 = strcat(root_f,cellstr(frame_names));
        fmriTemp =  spm_read_vols(spm_vol(frame_names3{1}));
        voxel_TCs = fmriTemp(rAICvoxels)'; 
        betaB = fmriTemp(rAICvoxels);%ROI voxels for the current beta image
        BB = nanmean(betaB);%average value over all voxels in the mask
        BaselineB(subject) = BB;
        clear fmriTemp
        
        %Upregulation regressor
        frame_namesU = spm_select('List', root_f, beta_images_B{2});%selects the upregulation regressor for baseline tun
        frame_names4 = strcat(root_f,cellstr(frame_namesU));
        fmriTemp =  spm_read_vols(spm_vol(frame_names4{1}));
        voxel_TCs = fmriTemp(rAICvoxels)'; 
        betaBU = fmriTemp(rAICvoxels);%ROI voxels for the current beta image
        BU = nanmean(betaBU);%average value over all voxels in the mask
        BaselineUp(subject) = BU;
        clear fmriTemp
        diffBaseline(subject,1) = BU-BB; %Regulation success during baseline run
        
        %%Transfer run
        %Baseline regressor
        root_f_transfer = [globalDir filesep 'sub-VP' subjDir filesep sessDir{2} filesep statsDir filesep];%Change to second NFB session
        frame_names_transB = spm_select('List', root_f_transfer, beta_images_T{1});
        frame_names_5 = strcat(root_f_transfer,cellstr(frame_names_transB));
        fmriTemp =  spm_read_vols(spm_vol(frame_names_5{1}));
        voxel_TCs = fmriTemp(rAICvoxels)'; 
        betaTB = fmriTemp(rAICvoxels);
        TB = nanmean(betaTB);
        TransferBaseline(subject) = TB;
        clear fmriTemp 
        
        %Upregulation regressor
        frame_names_transU = spm_select('List', root_f_transfer, beta_images_T{2});
        frame_names_6 = strcat(root_f_transfer,cellstr(frame_names_transU));
        fmriTemp =  spm_read_vols(spm_vol(frame_names_6{1}));
        voxel_TCs = fmriTemp(rAICvoxels)'; 
        betaTU = fmriTemp(rAICvoxels);%ROI voxels 
        TU = nanmean(betaTU);
        TransferUp(subject) = TU;
        diffTransf(subject,1) = TU-TB; %Regulation Success transfer run
        
        id{subject,1} =  {subjDir};
        diffSessions(subject,1) = diffTransf(subject,1)-diffBaseline(subject,1);%Overall NFB success per subject
  
        
end 

        
tableConn = table(id, diffBaseline,diffTransf,diffSessions)
writetable(tableConn,'ConnrAICBaseline_Transfer_fullDuration.xlsx') 
