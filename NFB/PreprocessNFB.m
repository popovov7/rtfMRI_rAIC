% Jeanette Popovova, last updated 07-08-2021
% general preprocessing pipeline for functional MRI images using SPM
% batches for Neurofeedback folder structure



clear all
clc

preprocFlag = 1;        % pre-process functional data
subj = [];       % enter subject list here

spm fmri                            % initiate spm
spm_specifics = what('spm12');      % get some spm info in a structure
SPMDir = spm_specifics.path;        % get the path to spm directory
addpath(SPMDir);                    % add path just in case
spm_jobman('initcfg')               % initiate the jobman

globalDir   = pwd;                 
rawFormat   = 'nii';         

%Specify directory names
statsDir = 'stats';
subjDir = 'CP_';%'VP' for experimental group
funcs = {'_ses-03_task-nfb_run-0_bold.nii' '_ses-03_task-nfb_run-1_bold.nii' '_ses-03_task-nfb_run-2_bold.nii' '_ses-03_task-nfb_run-3_bold.nii' '_ses-03_task-nfb_run-4_bold.nii' '_ses-03_task-nfb_run-5_bold.nii' '_ses-03_task-nfb_run-6_bold.nii'};%''_ses-02_task-nfb_run-0_bold.nii' '_ses-02_task-nfb_run-1_bold.nii' _ses-02_task-nfb_run-2_bold.nii' '_ses-02_task-nfb_run-3_bold.nii' '_ses-02_task-nfb_run-4_bold.nii' '_ses-02_task-nfb_run-5_bold.nii' '_ses-02_task-nfb_run-6_bold.nii' %names of the func files (without subject name)
struc = '_ses-03_T1w.nii';
onsetDir = 'onsets';
onset_prefix = 'onsets_r';  
run = {'r0' 'r1' 'r2' 'r3' 'r4' 'r5' 'r6'};
SessDir = 'ses-03';
runDir ='r';
structDir   = 'anat';%directory with structural image
funcDir = 'func';% directory with functional images
smoothK = [8 8 8];% desired smoothing kernel 
voxelRes = [1.97 2 3];
nr_slices = 27;             % how many slices in volume
TR = 2;                 % repition time (sampling frequency)
TA = TR - (TR/nr_slices);   % gap time between volumes
SO = [1:2:nr_slices 2:2:nr_slices];%slice order for Philipps scanner default interleaved ascending
RS = 1; 

% Advised order by UCL
% 1 – slice timing correction of EPI data       (not customary whe using multiband or quick TR (< 2 sec)
% 2 – realignment of EPI data                   (write mean only)
% 3 – coregister T1 to mean EPI
% 4 – segment the coregistered T1
% 5 – normalize EPI data (and possibly T1) using the deformation
%     field obtained by the segmentation
% 6 – smooth the functional data
% (7 – stats modelling and estimation)
% (8 - creating contrast images at the first level)


steps = {'slice_timing' 'realign' 'coreg' 'segment'  'normalize' 'smooth'}; 

for ss= subj
    subjDir = [sprintf('%03d',ss)];
    subj2Dir = [globalDir filesep 'sub-CP' subjDir]
    subjStructDir   = [subj2Dir filesep SessDir filesep structDir filesep];

        if preprocFlag == 1

            for ii = 1:length(steps)

                switch steps{ii}

                    case 'slice_timing'
                        
                        clear matlabbatch
                        cme = 1;
                       
                            
                            root_f = [subj2Dir filesep SessDir filesep funcDir filesep] % path to current subjects' func folder

                             for k=1:numel(run)
                                
                                root_f = [subj2Dir filesep SessDir filesep funcDir filesep run{k} filesep]                                                   
                                f1   = spm_select('List', root_f, ['sub-CP' subjDir funcs{k}]); 
                                f2  = cellstr([repmat(root_f,size(f1,1),1) f1]);                   
                                matlabbatch{1}.spm.temporal.st.scans(cme) = {f2};                        
                                if isempty(f2{1})   
                                    fprintf('\nno functional files loaded in: %s\n', steps{ii})     
                                end
                                cme=cme+1;
                            end
                       

                        % variables specified above
                        matlabbatch{1}.spm.temporal.st.nslices = nr_slices;
                        matlabbatch{1}.spm.temporal.st.tr = TR;
                        matlabbatch{1}.spm.temporal.st.ta = TA;
                        matlabbatch{1}.spm.temporal.st.so = SO;
                        matlabbatch{1}.spm.temporal.st.refslice = RS;
                        matlabbatch{1}.spm.temporal.st.prefix = 'a';        

                        % save the batch, run it and clear it
                        save([subj2Dir filesep 'slice_timing'], 'matlabbatch');
                        spm_jobman('run', matlabbatch);                     
                        clear matlabbatch
                        fprintf(['slice-timing DONE for subject: ' num2str(ss) '\n'])
                        
                    case 'realign'
                        
                        clear matlabbatch
                        cme=1;
                        root_f = [subj2Dir filesep SessDir filesep funcDir filesep];
                        
                            for k=1:numel(run)
                               
                                root_f = [subj2Dir filesep SessDir filesep funcDir filesep run{k} filesep];
                                f1   = spm_select('List', root_f, ['asub-CP' subjDir funcs{k}]); % will update for each run
                                f2  = cellstr([repmat(root_f,size(f1,1),1) f1]);
                                matlabbatch{1}.spm.spatial.realign.estwrite.data{cme} = f2; 
                                cme=cme+1;
                            end
                       
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

                        % save the batch, run it and clear it
                        save([subj2Dir filesep 'realign'], 'matlabbatch')
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        fprintf(['realignment (est/write) DONE for subject: ' num2str(ss) '\n'])

                    case 'coreg'

                        % retrieve the structural
                        f1struct        = spm_select('List', subjStructDir, ['sub-CP' subjDir struc]);
                        f2struct        = cellstr([repmat(subjStructDir,size(f1struct,1),1) f1struct]);

                        % retrieve the mean EPI image
                        root_f = [subj2Dir filesep SessDir filesep funcDir filesep run{1} filesep];%the mean image is saved in run 1, if multiple sessions udapt to session1
                        mean_image      = spm_select('List', root_f, ['^mean' '.*\.' rawFormat '$'])

                        % define reference (mean EPI) and source (struct)
                        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile([root_f filesep mean_image])};
                        matlabbatch{1}.spm.spatial.coreg.estimate.source = f2struct;
                        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
                        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

                        % save the batch, run it and clear it
                        save([subj2Dir filesep 'coregistration'], 'matlabbatch')
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        fprintf(['coregistration DONE for subject: ' num2str(ss) '\n'])

                    case 'segment'

                        % get the structural
                        f1struct        = spm_select('List', subjStructDir, ['sub-CP' subjDir struc]);
                        f2struct        = cellstr([repmat(subjStructDir,size(f1struct,1),1) f1struct]);
                        
                        matlabbatch{1}.spm.spatial.preproc.channel.vols(1) = f2struct;
                        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
                        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
                        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];

                        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,1']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
                        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,2']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
                        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,3']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
                        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,4']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
                        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,5']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
                        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[SPMDir filesep 'tpm' filesep 'TPM.nii,6']};
                        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
                        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
                        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];

                        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
                        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
                        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
                        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
                        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
                        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
                        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

                        % save the batch, run it and clear it
                        save([subj2Dir filesep 'segmentation'], 'matlabbatch')
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        fprintf(['Segmentation DONE for subject: ' num2str(ss) '\n'])


                    case 'normalize'

                        % ------- Normalize using the forward deformation file y_ from segmentation step
                        fd1        = spm_select('List', subjStructDir, ['y_sub-CP' subjDir struc]);
                        fd2        = cellstr([repmat(subjStructDir,size(fd1,1),1) fd1]);
                        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {};
                        root_f = [subj2Dir filesep SessDir filesep funcDir filesep];

                            for k=1:numel(run)
                                root_f = [subj2Dir filesep SessDir filesep funcDir filesep run{k} filesep];
                                f1   = spm_select('List', root_f, ['rasub-CP' subjDir funcs{k}])
                                f2  = cellstr([repmat(root_f,size(f1,1),1) f1]);
                                matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = fd2;
                                matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [matlabbatch{1}.spm.spatial.normalise.write.subj.resample; f2];
                            end
                      

                        save([subj2Dir filesep 'normalize'], 'matlabbatch')
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        fprintf(['Normalization DONE for subject: ' num2str(ss) '\n'])

                    case 'smooth'

                        matlabbatch{1}.spm.spatial.smooth.data = [];

                            for k=1:numel(run)  
                                root_f = [subj2Dir filesep SessDir filesep funcDir filesep run{k} filesep]
                                f1   = spm_select('List', root_f, '.*^wra')
                                f2  = cellstr([repmat(root_f,size(f1,1),1) f1]);
                                matlabbatch{1}.spm.spatial.smooth.data = [matlabbatch{1}.spm.spatial.smooth.data; f2];
                            end

                        matlabbatch{1}.spm.spatial.smooth.fwhm = smoothK;
                        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                        matlabbatch{1}.spm.spatial.smooth.im = 0;
                        matlabbatch{1}.spm.spatial.smooth.prefix = 's2';

                        save([subj2Dir filesep 'smooth'], 'matlabbatch')
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch
                        fprintf(['Smoothing DONE for subject: ' num2str(ss) '\n'])
                end
       


            end
        end

end
