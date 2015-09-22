%% getUnits
% The script gets all chirp units with a specified ranksum (Miro's correlation criterion)

% Find all chirp keys for LGN with sorted units
chirp_keys = fetch(data.Series & data.Experiments('exp_name LIKE "%chirp%"'));

% Exclude all experiments before '2014-03-14'
date_keys = fetch(data.Series(chirp_keys),'series_date');
mask_date = datenum(fetchn(data.Series(chirp_keys),'series_date')) > datenum('2014-03-14', 'yyyy-mm-dd');

% Contains all valid chirp keys
chirp_keys = date_keys(mask_date);
chirp_keys = rmfield(chirp_keys,'series_date');

% corr_p represents the wilcoxon ransum correlation for between/withing
% segments
units_for_chirp = fetch(miro.ChirpQuality(chirp_keys) & miro.ChirpQuality('corr_p < 0.00001'));
