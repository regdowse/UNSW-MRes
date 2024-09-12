%MAIN_AMEDA
%
%   MAIN_AMEDA is the main function of the eddy detection and
%   tracking package. It returns position of the centers, dimensions and 
%   tracks of the eddies detected from the time series of a 2-D velocity 
%   field.
%   It gives also an history of the splitting and merging events.
%
%   - 'source' allows to specify the type of sources file (AVISO, ROMS, NEMO,...)
%     with their specific parameters and Input/Output.
%   - cpus to use 'parfor' as time loops (=# of processors)
%       cpus = 1 (default)
%   - update is a flag allowing to update an existing tracking:
%       update = number of time steps backward to consider
%       update = 0 (default) to compute all the time serie
%   - stepF is the last time step computed
%       stepF = temporal size of the input data
%
%   The algortihm subroutines:
%
%   - mod_eddy_params sets user defined paths and parameters:
%     nc_u nc_v nc_dim b bx r path_in path_out periodic criteres
%     Users should modify keys_sources.m according to their 
%     settings.
%
%   - mod_init initialise or update mat-file.
%
%   - mod_fields compute LNAM.
%
%   - mod_eddy_centers returns a structure array with the position of the
%     detected eddy centers.
%
%   - mod_eddy_shapes computes dimensions for the detected eddy centers.
%
%   - mod_eddy_tracks computes eddy tracks using the detected centers.
%
%   Find the output files in path_out:
%
%   - fields.mat contains detection_fields with LNAM for each step.
%   - eddy_centers.mat contains for each step:
%       * centers0 as the local max(LNAM)
%       * centers as the potential centers
%       * centers2 as the detected eddies
%   - eddy_shapes.mat contains for each step:
%       * shapes1 the eddy features
%       * shapes2 the common double contour features
%       * profil2 the streamlines features scanned around each eddy
%       * warn_shapes the flag for potential centers
%       * warn_shapes2 the flag for detected eddies
%   - eddy_tracks.mat contains eddy centers, features and flags for each eddy
%
%-------------------------
%   June 2016 Briac Le Vu
%-------------------------
%
%=========================

start
clear; clc;

%----------------------------------------
% source of data driving the netcdf format
source = 'AVISO';

%----------------------------------------
% domaine
keys = 'MED';

%----------------------------------------
% Update option
update = 0; % the serie from the begenning

%----------------------------------------
% Possibility to shorter the serie
%stepF = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Produce default parameters in param_eddy_tracking
if exist('stepF','var')
    mod_eddy_params(['keys_sources_',source,'_',keys],stepF)
else
    mod_eddy_params(['keys_sources_',source,'_',keys])
end
run(['keys_sources_',source,'_',keys])
load('param_eddy_tracking','path_out','streamlines','resol','stepF')

%----------------------------------------
% Preallocate structure array and mat-file or prepare update
% !! replace or reinitialise previous results !!
step0 = mod_init(stepF,update);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main routines ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Build I/O matfile
disp('Your MATLAB to support "-v7.3" format to get full performance of')
disp('I/O MAT-file and save memory space')
disp(' ')

%----------------------------------------
% Build I/O matfile
fields_inter_mat = matfile([path_out,'fields_inter.mat'],'Writable',true);
fields_mat = matfile([path_out,'fields.mat'],'Writable',true);
centers_mat = matfile([path_out,'eddy_centers.mat'],'Writable',true);
shapes_mat = matfile([path_out,'eddy_shapes.mat'],'Writable',true);

%----------------------------------------
% prepare log directory
rmdir([path_out,'/log'])
mkdir([path_out,'/log'])

for stp = step0:stepF
    
    %----------------------------------------
    % begin the log file
    if stp<10
        diary([path_out,'log/log_eddy_stp_000',num2str(stp),'.txt']);
    elseif stp<100
        diary([path_out,'log/log_eddy_stp_00',num2str(stp),'.txt']);
    elseif stp<1000
        diary([path_out,'log/log_eddy_stp_0',num2str(stp),'.txt']);
    elseif stp<10000
        diary([path_out,'log/log_eddy_stp_',num2str(stp),'.txt']);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute LNAM ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Compute non interpolated fields for step stp
    fields = mod_fields(source,stp,1);
    
    %----------------------------------------
    % Write in I/O matfile step stp
    fields_mat.detection_fields(1,stp) = fields;
    
    if resol>1
        %----------------------------------------
        % Compute interpolated fields for step stp
        fields = mod_fields(source,stp,resol);
    end
    
    %----------------------------------------
    % Write in I/O matfile
    % Interpolated and non interpolated field can be the same
    fields_inter_mat.detection_fields(1,stp) = fields;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find centers ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Detection of eddy centers for step stp
    [centers0,centers] = mod_eddy_centers(source,stp,fields);
    
    %----------------------------------------
    % Write in I/O matfile step stp
    centers_mat.centers0(1,stp) = centers0;
    centers_mat.centers(1,stp) = centers;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------------------------
    % Determination of eddy features for step stp
    if streamlines
        [centers2,shapes1,shapes2,profil2,warn_shapes,warn_shapes2] = ...
            mod_eddy_shapes(source,stp,fields,centers);
    else
        [centers2,shapes1,shapes2,~,warn_shapes,warn_shapes2] = ...
            mod_eddy_shapes(source,stp,fields,centers);
    end
    
    %----------------------------------------
    % Write in I/O matfile step stp
    centers_mat.centers2(1,stp) = centers2;
    shapes_mat.shapes1(1,stp) = shapes1;
    shapes_mat.shapes2(1,stp) = shapes2;
    shapes_mat.warn_shapes(1,stp) = warn_shapes;
    shapes_mat.warn_shapes2(1,stp) = warn_shapes2;
    if streamlines
        shapes_mat.profil2(1,stp) = profil2;
    end
    
    %----------------------------------------
    % close log file
    diary off

end

%----------------------------------------
% concatene log files
if exist([path_out,'log_eddy_stp.txt'],'file')
    movefile([path_out,'log_eddy_stp.txt']);
end
system(['cat ',path_out,'log/log_eddy_stp*.txt >> ',path_out,'log_eddy.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track eddies ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% Tracking eddies and record interacting events
mod_eddy_tracks_nopool('',update)

%----------------------------------------
% Resolve merging and spltting event and filter eddies shorter than cut_off
mod_merging_splitting('');

