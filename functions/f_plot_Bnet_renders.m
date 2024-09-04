function nodal_strength = f_plot_Bnet_renders(T_mask,path2Results,FileName,flag_mean,Path2Parc)    
BrainNet_scripts = fullfile('..','functions/BrainNetViewer'); % For nodal strenght nifti files to open in BrainNet 
addpath(genpath(BrainNet_scripts));
MRI_scripts = fullfile('..', 'functions/','toolbox_matlab_nifti'); % Needed to load nifti in matlab
addpath(genpath(MRI_scripts));
Path2fig = fullfile(path2Results,'brain_renders');
GM_parc = load_nifti(fullfile(Path2Parc,'shen_dil.nii.gz'));                       
if ~isfolder(Path2fig)
    mkdir(Path2fig);
end

GM_parc_temp = GM_parc;
nReg = max(GM_parc_temp.vol(:));       
if flag_mean == 1
    nodal_strength = nanmean(T_mask);
elseif flag_mean == 2
    nodal_strength = nansum(T_mask);
elseif flag_mean == 3
    nodal_strength = nansum(abs(T_mask));
else
    nodal_strength = T_mask;
end
RC_surface = zeros(size(GM_parc_temp.vol));            
for j=1:nReg
    RC_surface(GM_parc_temp.vol==j) = nodal_strength(j);
end
GM_parc_temp.vol = RC_surface;  
save_nifti(GM_parc_temp,fullfile(Path2fig,[FileName 'shen_dil.nii.gz']));
