% Jeanette Popovova, last updated 28-08-2020
% This script will specify, estimate and build the contrast for NFB
% sessions


clear all
clc
spm fmri                            % initiate spm
spm_specifics = what('spm12');      % get some spm info in a structure
SPMDir = spm_specifics.path;        % get the path to spm directory
addpath(SPMDir);                    % add path just in case
spm_jobman('initcfg')               % initiate the jobman

stats_1stlevel = 1;   % performing stats at the first level putting several sessions in the model (if you have more than one)
estimation = 1;    % do want to specify and etimate the model?    
contrast  = 1;    % do you want to specify and create the conof interest images?

subj = [];       % enter subject list here
globalDir   = pwd;                 
rawFormat   = 'nii';                
statsDir = 'statsNewJun';
subjDir = 'VP';
funcsName = {'_task-nfb_run-0_bold.nii' '_task-nfb_run-1_bold.nii' '_task-nfb_run-2_bold.nii' '_task-nfb_run-3_bold.nii' '_task-nfb_run-4_bold.nii' '_task-nfb_run-5_bold.nii' '_task-nfb_run-6_bold.nii'}
onset_name = 'OnsetsNFBSep';  %name of Onset file
runDir = {'r0' 'r1' 'r2' 'r3' 'r4' 'r5' 'r6'};%the runs have to be saved in seperate folders
sessDir = {'ses-02'}; %specify which session
structDir   = 'anat';%directory with structural image
funcDir = 'func';% directory with functional images inside runfolder
TR=2;


for ss= subj
    subjDir = [sprintf('%03d',ss)];
    subj2Dir = [globalDir filesep 'sub-VP' subjDir]
    subjStructDir   = [subj2Dir filesep sessDir filesep structDir filesep];

    for session=1:numel(sessDir)
        Modeldir =  [subj2Dir filesep sessDir{session} filesep statsDir];
        
        if stats_1stlevel == 1
        counter =1; 
        matlabbatch{1}.spm.stats.fmri_spec.dir = {Modeldir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 33;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        for k=1:numel(runDir)
            funcDir     = [subj2Dir filesep sessDir{session} filesep 'func' filesep runDir{k} filesep];     
            onsetsDir =  [globalDir filesep onset_name '.mat'];
            rp1 = spm_select ('List', funcDir, ['^rp']);%motion regressors
            rp2 = cellstr([repmat(funcDir,size(rp1,1),1) rp1]);
            funcsName2= ['s2wrasub-VP' subjDir '_' sessDir{session} funcsName{k}]
            f1= spm_select('expand', fullfile(funcDir,funcsName2))
            %f1   = spm_select('List', funcDir, ['^ss' '.*\.' rawFormat '$']);
            f2  = cellstr([f1]);  
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).scans = f2;
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).multi = {onsetsDir};
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).multi_reg = rp2;
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).hpf = 128;

            counter=counter+1;
             
         end
            
         
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
          
            spm_jobman('run', matlabbatch);
            clear matlabbatch
         
            fprintf(['modeling DONE for subject: ' num2str(ss) '\n'])
          
             end
            
        
        if estimation==1

             matlabbatch{1}.spm.stats.fmri_est.spmmat = {[Modeldir filesep 'SPM.mat']};
             matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
             matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

             spm_jobman('run', matlabbatch);
             clear matlabbatch

                fprintf(['modeling DONE for subject: ' num2str(ss) '\n'])   
         end
    
        if contrast == 1
            clear matlabbatch 
            statsdir =  [subj2Dir filesep sessDir{session} filesep statsDir];          
            matlabbatch{1}.spm.stats.con.spmmat = {[statsdir filesep 'SPM.mat']};

            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'r0';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0 zeros(1,6)];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'r1';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights =  [zeros(1,9) -1 1 0 zeros(1,6)];
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'r2';
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [zeros(1,18) -1 1 0 zeros(1,6)];
            matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'r3';
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [zeros(1,27) -1 1 0];
            matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'r4';
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [zeros(1,36) -1 1 0];
            matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'r5';
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [zeros(1,45) -1 1 0];
            matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'r6';
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,54) -1 1 0];
            matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;

            matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'All';
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = repmat([-1 1 0 zeros(1,6)],1,7);
            matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;
            spm_jobman('run',matlabbatch); 
            clear matlabbatch 

        end
         end
        end

