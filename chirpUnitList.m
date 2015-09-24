function chirpUnitList = chirpUnitList()
%% chirpUnitList = chirpUnitList()
% The script creates a table of all units that satisfy specific quality criteria in their response
% to the chirp stimulus.

%% Startup datajoint
startup_cin

%% Parameters
corr_p = 0.0001; % corr_p represents the wilcoxon ransum correlation for between/withing segments
qi = 0.075; % Philipps quality index
miro_qi = 3; % Miros quality index of likely single- vs multiunit

%% Find all chirp keys for LGN with sorted units
chirp_keys = fetch(data.Series & data.Experiments('exp_name LIKE "%chirp%"'));

% Exclude all experiments before '2014-03-14'
date_keys = fetch(data.Series(chirp_keys),'series_date');
mask_date = datenum(fetchn(data.Series(chirp_keys),'series_date')) > datenum('2014-03-14', 'yyyy-mm-dd');

% Contains all valid chirp keys
chirp_keys = date_keys(mask_date);
chirp_keys = rmfield(chirp_keys,'series_date');

% Create table of units that respond satisfactorily to chirp stimulus
chirpUnitList = fetch(miro.ChirpQuality(chirp_keys) ...
    & miro.ChirpQuality(sprintf('corr_p < %f',corr_p)) ...
    & miro.ChirpQuality(sprintf('berens_qi >= %f',qi)) ...
    & data.ClusterInfo(sprintf('quality <= %f',miro_qi)),...
    'corr_p', 'berens_qi');


end