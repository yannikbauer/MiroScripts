%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LGN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FLAGS & VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up
%clear all;
%startup;

% Visual Sorting
flag.contrast  = 0;    % 1 = enable; 0 = disable
flag.ranksum   = 1;    % 1 = enable; 0 = disable

spkNumThreshold  = 300;  % less than # will be excluded from the analysis

% Stim flag
flag.chirp     = 1;    % Chirp
flag.ac        = 0;    % AcTun
flag.sparse    = 0;    % SparseNoise4Blank
flag.aori      = 0;    % AoriTun
flag.atf       = 0;    % AtfTun
flag.size      = 0;    % SizeTuning
flag.ds        = 0;    % DS Bar
flag.asw       = 0;    % Average Spike Wave
flag.autoc     = 0;    % Autocorrelogram

% Results
flag.kill      = 1;    % 1 = kills all data in the results folder
flag.save      = 1;    % 1 = saves results plots and CLOSES ALL FIGURES
flag.populate  = 0;    % 1 = activates automatic populating of empty files


% Override Unit
plot_iunit     = 0;    % 0 = disable, i.e. plots all units

% Choose mouse name (BL6_0195), if you want to process a mouse separately
mouse_id       = '';      % '' for all mice

% Pathway for saving the results
pathway = '/Users/miro/Desktop/Figures_xx';
if(exist(pathway, 'dir') ~= 7)
    mkdir(pathway)
end

% Data points need to accept a unit as visually driven
nVis = 2;





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INITIALIZE STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear results will delete all files in the results folder
if(flag.kill == 1)
    rmdir(pathway, 's');
    mkdir(pathway)
end

% Open Txt File for writing error messages
stim_exceptions = fopen(fullfile(pathway,'_exceptions'), 'wt');
fprintf(stim_exceptions,'%s\t%s\t%s\t%s\t%s\t%s\n', 'Stimulus', 'Message', 'mouse_counter', 'series_num','exp_num','unit_id');






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  KEY DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get chirp keys for one single mouse or all
if(isempty(mouse_id) == 0)
    chirp_keys.mouse_counter = fetch1(data.Mice(sprintf('mouse_id LIKE "%s"', mouse_id)), 'mouse_counter');
    chirp_keys = fetch(data.Experiments(chirp_keys) & 'exp_name LIKE "%Chirp%"');
else
    % Find all chirp keys for LGN with sorted units
    chirp_keys = fetch(data.Experiments('exp_name LIKE "%chirp%"'));
end

% Get the series keys
skeys = fetch(data.Series & data.Experiments(chirp_keys));





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  POPULATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It checks again, whether everything is populated
if(flag.populate)
    counter = 1;  % exception count
    for ikey = 1 : numel(chirp_keys)
        try
            populate(data.Stimuli, chirp_keys(ikey));
            populate(data.Channels, chirp_keys(ikey));
            populate(data.ClusterInfo, chirp_keys(ikey));
            populate(data.TrialSpikes, chirp_keys(ikey));
            populate(data.TrialSpikesExtra, chirp_keys(ikey));
            populate(data.ConditionSpikes, chirp_keys(ikey));
            populate(data.ConditionSpikesExtra, chirp_keys(ikey));
            populate(data.StimInfo, chirp_keys(ikey));
            populate(miro.ChirpQuality, chirp_keys(ikey));
            populate(data.Tuning, chirp_keys(ikey));
            populate(data.Locomotion, chirp_keys(ikey));
            populate(data.SeriesEvents, chirp_keys(ikey));
            populate(data.TrialSDFs, chirp_keys(ikey));
            populate(data.ConditionSDFs, chirp_keys(ikey));
            populate(data.SpikeQuality, chirp_keys(ikey));
        catch ME
            display('skipped - populate error');
            mouse_id = fetch1(data.Mice(chirp_keys(ikey)), 'mouse_id');
            
            % Save error in the txt-file in the figure folder
            fprintf(stim_exceptions, '%s\t', 'Chirp');
            fprintf(stim_exceptions, '%s\t', 'Not populated');
            fprintf(stim_exceptions, '%s\t', mouse_id);
            fprintf(stim_exceptions, '%d\t', chirp_keys(ikey).series_num);
            fprintf(stim_exceptions, '%d\t', chirp_keys(ikey).exp_num);
            fprintf(stim_exceptions, '%d\t\n', -1);
            
            % save exceptions as mat file
            MExceptions(counter).key.mouse_counter = chirp_keys(ikey).mouse_counter;
            MExceptions(counter).key.series_num = chirp_keys(ikey).series_num;
            MExceptions(counter).key.exp_num = chirp_keys(ikey).exp_num;
            MExceptions(counter).identifier = ME.identifier;
            MExceptions(counter).message = ME.message;
            MExceptions(counter).cause = ME.cause;
            MExceptions(counter).stack = ME.stack;
            
            counter = counter + 1;
            
            % Check for override
            if(plot_iunit > 0)
                error('Populating failed');
            else
                continue;
            end
        end
    end
    
    % If any exceptions, save them
    if(exist('MExceptions','var') && isstruct(MExceptions))
        path = fullfile(pathway, 'Data');
        if(exist(path, 'dir') ~= 7)
            mkdir(path)
        end
        save(fullfile(path, 'MExceptions.mat'), 'MExceptions');
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  UNITS FOR CHIRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get units for chirps
units_for_chirp = fetch(data.Units & data.Experiments(chirp_keys));
[units_for_chirp.exp_num_chirp] = units_for_chirp.exp_num;
units_for_chirp = rmfield(units_for_chirp,'exp_num');

% Initialize contrast, ranksum & qi fields
[units_for_chirp.contrast] = deal(-1);
[units_for_chirp.ranksum] = deal(100);
[units_for_chirp.qi] = deal(0);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1. DETERMINATION OF VISUAL RESPONSIVENESS - CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The responsiveness of each unit for the chirp stimulus will be tested
% with the responsiveness for the contrast tuning stimulus (AcTun)

% Disable or Enable
if(flag.contrast == 1)
    
    % Get AcTun & Chirp keys
    vkeys = fetch(data.Experiments(skeys) & 'exp_name LIKE "%AcTun%"');
    
    % Initialize
    contrast_key = [];
    
    % Go through all visual keys and determine whether the units are
    % responsive to chirp and AcTun stimuli
    for iseries = 1 : numel(vkeys)
        fprintf('Computing contrast %d of %d\n', iseries, numel(vkeys));
        ukeys = fetch(data.Units(vkeys(iseries)), 'unit_id'); % contains exp_num & units
        stimInfo = fetch(data.StimInfo(vkeys(iseries)), '*');
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Visual responsivenss %
        %%%%%%%%%%%%%%%%%%%%%%%%
        uLog = false(1, numel(ukeys));
        for iunit = 1 : numel(ukeys)
            c = fetch(data.ConditionSpikes(ukeys(iunit)) & data.GratingConditions('grat_opto_light = 0'), 'grat_num', 'cond_rate', 'cond_sem');
            activeLog = ismember([c.grat_num], stimInfo.num_active_grats);
            blankLog  = ismember([c.grat_num], stimInfo.num_blank_grats);
            % positive response || negative response
            if nnz([c(activeLog).cond_rate] - 2.58*[c(activeLog).cond_sem] > c(blankLog).cond_rate) > nVis || ...
                    nnz([c(activeLog).cond_rate] + 2.58*[c(activeLog).cond_sem] < c(blankLog).cond_rate) > nVis
                uLog(iunit) = true;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign experiment number %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(uLog)
            % remove experiment number
            temp = rmfield(ukeys(uLog), 'exp_num');
            % find the relevant chirp experiment
            chirp_log = [chirp_keys.mouse_counter] == temp(1).mouse_counter & [chirp_keys.series_num] == temp(1).series_num;
            
            
            %%%%%%%%%%%%%%%%%%%
            % MULTIPLE CHIRPS %
            %%%%%%%%%%%%%%%%%%%
            % In the case we have more than one chirp stimulus/series,
            % compare the monitor angle with AcTun8 otherwise take the first one
            if nnz(chirp_log)>1
                
                % Find the relevant experiments
                ckey_log = [vkeys.mouse_counter] == temp(1).mouse_counter & [vkeys.series_num] == temp(1).series_num;
                temp_ac = vkeys(ckey_log);          % keys for contrast experiment
                temp_chirp = chirp_keys(chirp_log); % keys for chirp experiment
                
                % Update temp_heys for monitor elevation and angle
                temp_ac = fetch(data.Experiments(temp_ac), 'exp_monitorangle', 'exp_monitorelevation');
                temp_chirp = fetch(data.Experiments(temp_chirp), 'exp_monitorangle', 'exp_monitorelevation');
                
                flag.monitor = 0;
                
                for i = 1:numel(temp_chirp) % all chirp keys
                    if (temp_chirp(i).exp_monitorangle == temp_ac.exp_monitorangle)  % same monitor angle
                        [temp.exp_num_chirp] = deal(temp_chirp(i).exp_num);
                        flag.monitor = 1;
                        break
                    end
                end
                
                % Different monitor angle in all cases, take the first one
                if(flag.monitor == 0)
                    temp = temp_chirp(1).exp_num;
                end
            else
                % There is only one chirp stimulus/experiment
                [temp.exp_num_chirp] = deal(chirp_keys(chirp_log).exp_num);
            end
            contrast_key = [contrast_key; temp]; %#ok<AGROW>
        end
    end
    
    % Extra variable; units_for_chirp and contrast have to have the same
    % number of fields
    contrast = zeros(numel(units_for_chirp),1);
    
    % Check, which structures passed the contrast test
    for icunit = 1:numel(contrast_key)
        for iunit = 1:numel(units_for_chirp)
            if(isequal(contrast_key(icunit), units_for_chirp(iunit)))
                contrast(iunit) = 1;
            end
        end
    end
    
    % Assign the contrast flag
    contrast = num2cell(contrast);
    [units_for_chirp.contrast]=deal(contrast{:});
end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2. DETERMINATION OF VISUAL RESPONSIVENESS - CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To determine whether a neuron was visually driven, we separated the movie
% into two segments (e.g., separated a 6-s movie into two 3-s segments),
% and computed the average between-trial CCs (responses binned at the
% stimulus frame rate) within segment and between segments. Only those
% cells that had significantly higher within-segment CCs (P < 0.01,
% Wilcoxon rank sum test) and firing rate > 1 spike/s were included in the
% analyses.

for iunit = 1 : numel(units_for_chirp)
    try
        % Plot all or one unit; 0 = disable
        if(plot_iunit > 0)
            iunit = plot_iunit;
        end
        
        fprintf('Computing ranksum for unit %d of %d\n', iunit, numel(units_for_chirp));
        
        % Define key
        thisUnit = units_for_chirp(iunit);
        thisUnit.exp_num = units_for_chirp(iunit).exp_num_chirp;
        
        % Get all trials
        [trials] = fetch(data.TrialSpikesExtra(thisUnit), 'trial_num');
        
        % get all trial onsets in samples
        [trial_onsets, trial_offsets] = arrayfun(@(key) fetch1(data.ExtraTrials(key), 'trial_onset', 'trial_offset'), trials);
        
        if(isempty(trial_onsets))
            errnum = 1;
            error('Trial_onsets error')
        end
        
        % get trial duration
        trial_durs = trial_offsets - trial_onsets;
        trial_duration  = median(trial_durs);
        seg1_start = 0;
        seg1_end = trial_duration/2;
        seg2_start= trial_duration/2 + 1;
        seg2_end = trial_duration;
        
        % get unit spiketimes
        unitSpkTs = fetch1(data.Units(thisUnit), 'unit_spiketimes');
        
        % low firing rate, continue - will be sorted out below
        if(size(unitSpkTs,1) < spkNumThreshold)
            display('skipped - low firing rate')
            if(plot_iunit > 0)
                break
            else
                continue
            end
        end
        
        % for each trial, we want to get the spike counts during the first and
        % second segments for each bin
        binSize = 30000/60; % binnig at monitor stimulus frame rate for 30000 samples per s and 60 frames per s
        
        % example trial
        nBinsPerSeg = floor((trial_duration/2)/binSize); % to get rid of the last bin in the histc
        S = struct('type', '()', 'subs', {{1:nBinsPerSeg}});  % for subsref
        
        % lets put it inside array fun
        seg1Trains = arrayfun(@(on, off) subsref(histc(unitSpkTs, on+seg1_start:binSize:on+seg1_end),S) , trial_onsets, trial_offsets, 'uniformoutput', 0) ;
        seg2Trains = arrayfun(@(on, off) subsref(histc(unitSpkTs, on+seg2_start:binSize:on+seg2_end),S) , trial_onsets, trial_offsets, 'uniformoutput', 0) ;
        
        % corr correlates between columns of a matrix (or matrices)
        % for within seg corr(seg1Trains) and corr(seg2Trains)
        % for between corr(seg1, seg2).
        seg1_corr = corr(cell2mat(seg1Trains'));
        seg2_corr = corr(cell2mat(seg2Trains'));
        seg1_seg2_corr = corr(cell2mat(seg1Trains'),cell2mat(seg2Trains'));  % diagonale nehmen
        
        % Withing:
        % Apply a logical mask for withing segment correlations where we
        % discard correlations of the same trials and repetition. The
        % values will be stored in one vector
        
        mask_corr = triu(ones(size(seg1_corr)),1);
        %seg_vec = cat(1, seg1_corr(mask_corr>0), seg2_corr(mask_corr>0)); 
        
        if mean(mean(seg1_corr(mask_corr>0))) > mean(mean(seg2_corr(mask_corr>0)))
            seg_vec = seg1_corr(mask_corr>0);
        else
            seg_vec = seg2_corr(mask_corr>0);
        end
       
        % Between:
        % in order to avoid general drift between separate trials in
        % between segment condition, take only coresponding trial-to-trial
        % values, i.e. the diagonal
        seg1_seg2_vec = diag(seg1_seg2_corr);
        
        % No spikes -> Test for NaNs
        if(isnan(nanmean(seg_vec)) == 1 || isnan(nanmean(seg1_seg2_vec)) == 1)
            errnum = 2;
            error('Correlation failed, contains only NaNs!');
        end
        
        % compute the Wilcoxon rank sum test for equal medians. H=1 indicates that
        % the null hypothesis can be rejected at the 5% level
        % [P,H] = ranksum(seg_vec, seg1_seg2_vec);
        units_for_chirp(iunit).ranksum = ranksum(seg_vec, seg1_seg2_vec);
                
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  3. DETERMINATION OF VISUAL RESPONSIVENESS - BERENS QUALITY CHECK
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Define key
        thisUnit = units_for_chirp(iunit);
        thisUnit.exp_num = units_for_chirp(iunit).exp_num_chirp;
        
        % to get rid of the last bin in the histc; here we have bins for the whole stimulus
        nBinsPerSeg = floor((trial_duration)/binSize);
        S = struct('type', '()', 'subs', {{1:nBinsPerSeg}});
        
        % lets put it inside array fun
        chirpTrains = arrayfun(@(on, off) subsref(histc(unitSpkTs, on:binSize:off),S) , trial_onsets, trial_offsets, 'uniformoutput', 0) ;
        
        % Philipps criterion for ca2+ imaging in RGCs; The threshold will have to
        % adjusted accordingly to these datafiles.
        f = @(d) var(mean(d,2),[],1)/mean(var(d,[],1),2);
        units_for_chirp(iunit).qi = f(cell2mat(chirpTrains'));
        
    catch
        
        display('skipped - exception handling');
        mouse_id = fetch1(data.Mice(thisUnit), 'mouse_id');
        
        % Save error in the txt-file in the figure folder
        fprintf(stim_exceptions, '%s\t', 'Chirp');
        if(errnum == 1)
            fprintf(stim_exceptions, '%s\t', 'Trial_onsets empty');
        elseif(errnum == 2)
            fprintf(stim_exceptions, '%s\t', 'NaNs in Corr');
        end
        fprintf(stim_exceptions, '%s\t', mouse_id);
        fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
        fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
        fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
        
        % Check for override
        if(plot_iunit > 0)
            if(errnum == 1)
                error('Trial_onsets empty');
            elseif(errnum == 2)
                errro('NaNs in Corr');
            end
        else
            continue;
        end
    end
    
    % Plot all or one unit; 0 = disable
    if(plot_iunit > 0)
        break
    end
end






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT CHIRP UNITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units for ranksum >=.05 and contrast == 0 will be rejected. Empty fields
% will be assign values so they could be discarded

if(plot_iunit == 0)
    mask_chirp = zeros(numel(units_for_chirp),1);
else
    mask_chirp = zeros(1,1);
end

if(flag.contrast || flag.ranksum)
    for iunit = 1:numel(units_for_chirp)
        
        % Plot all or one unit; 0 = disable
        if(plot_iunit > 0)
            iunit = plot_iunit;
        end
        
        fprintf('Sorting chirp units: %d of %d\n', iunit, numel(units_for_chirp));
        
        % Define key
        thisUnit = units_for_chirp(iunit);
        thisUnit.exp_num = units_for_chirp(iunit).exp_num_chirp;
        
        %%%%%%%%%%
        % VISUAL %
        %%%%%%%%%%
        % Sort for contrast or ranksum
        if(thisUnit.ranksum < 0.05 && flag.ranksum) || (thisUnit.contrast == 1 && flag.contrast)
            mask_chirp(iunit) = 1;
        end
        
        %%%%%%%%%%%%%%%
        % FIRING RATE %
        %%%%%%%%%%%%%%%
        % get unit spiketimes
        unitSpkTs = fetch1(data.Units(thisUnit), 'unit_spiketimes');
        
        % Get rid of less spiking units. If the unit passed the previous
        % test, then set the mask to zero again
        if(size(unitSpkTs,1) < spkNumThreshold)
            mask_chirp(iunit) = 0;
        end
        
        % Plot all or one unit; 0 = disable
        if(plot_iunit > 0)
            break
        end
    end
    
    % Apply the sorting mask
    units_for_chirp_sorted = units_for_chirp(mask_chirp == 1);
else
    
    % DO NOT SORT
    units_for_chirp_sorted = units_for_chirp;
end






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN CHIRP UNITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of chirp experiments with wrong stimulus, i.e. wrong monitor
% frequency, MID1s or trigger

mask_chirp_ws = ones(numel(units_for_chirp_sorted),1);

for iunit = 1:numel(units_for_chirp_sorted)
    
    fprintf('Cleaning chirp units: %d of %d\n', iunit, numel(units_for_chirp_sorted));
    
    % Define key
    thisUnit = units_for_chirp_sorted(iunit);
    thisUnit.exp_num = thisUnit.exp_num_chirp;
    
    %%%%%%%%%%
    % MID 1s %
    %%%%%%%%%%
    % Get rid of experiments with the wrong MID 1s stim part 
    series_date = fetch1(data.Series(thisUnit), 'series_date');
    if datenum(series_date, 'yyyy-mm-dd') <= datenum('2014-03-14', 'yyyy-mm-dd')
        mask_chirp_ws(iunit) = 0;
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WRONG MONITOR FREQUENCY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % expo data set for the stimulus
    mouse_id = fetch1(data.Mice(thisUnit), 'mouse_id');
    [~,~,unsortedFolder] = getEphysDataFolder(mouse_id, thisUnit.series_num);
    xmlFile = dir(fullfile(unsortedFolder, sprintf('*_%02d.xml', thisUnit.exp_num)));
    
    if(isempty(unsortedFolder))
        mask_chirp_ws(iunit) = 0;
        continue
    end
    
    % get the monitor frequency
    [~, expo_line] = unix(sprintf('grep TickDuration= %s', fullfile(unsortedFolder, xmlFile.name)));
    frameDur = str2double(regexp(expo_line, '(?<=TickDuration=")\d+', 'match'))/1000;
    
    
    if isempty(frameDur)
        continue
    end
    
    
    
    
    monitor_refreshRate = round(1/frameDur*1000);
    if(monitor_refreshRate == 120)
        mask_chirp_ws(iunit) = 0;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WRONG MONITOR FREQUENCY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % up to BL6_0276
    
end

% Apply the sorting mask
units_for_chirp_presorted = units_for_chirp_sorted;
units_for_chirp_sorted = units_for_chirp_sorted(mask_chirp_ws == 1);






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD OTHER STIMULI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stimuli definition
sStimList = {'cTun' 'oriTun' 'tfTun' 'sizeTun' 'sparseNoise4Blank' 'sfTun' 'DS'};
[~,nStimSize] = size(sStimList);

for istim = 1:nStimSize
    
    % Create stim string
    fprintf('Adding exp_num for %s\n', sStimList{istim});
    exp_num_var = sprintf('exp_num_%s',sStimList{istim});
    
    % Add exp_num to the key
    for iunit = 1:numel(units_for_chirp_sorted)
        units_for_chirp_sorted(iunit).(exp_num_var) = fetchn(data.Units(units_for_chirp_sorted(iunit)) & data.Experiments(sprintf('exp_name LIKE "%%%s%%" AND exp_name NOT LIKE "%%opto%%"', char(sStimList(istim)))), 'exp_num');
    end
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check, whether there is something to be plotted
if numel(units_for_chirp_sorted) == 0
    error('analyzeChirp:Plot','Units_for_chirp contains no units, plotting aborted')
end

% Layout Properties
[figPars, axPars] = setPlotPars('slide');
axPars.Units = 'centimeter';

% Initialize var
chirp_sdfs = cell(1,numel(units_for_chirp_sorted));

% Start plotting units
for iunit = 1 : 10%numel(units_for_chirp_sorted)
    
    
    
    
    % Plot all or one unit; 0 = disable
    if(plot_iunit > 0 && numel(units_for_chirp_sorted) > 1)
        iunit = plot_iunit;
    end
    
    close all
    
    % Open figure
    h = 42;   % height of the figure, specify for figures
    fh = figure(figPars, 'Position', [10 10 29.7 h], 'PaperUnits', 'centimeters', 'PaperType', 'a3');
    fprintf('Plotting unit %d of %d\n', iunit, numel(units_for_chirp_sorted));
    
    % Define key
    thisUnit = units_for_chirp_sorted(iunit);
    mouse_id = fetch1(data.Mice(thisUnit), 'mouse_id');
    latex_mouse_id = strrep(mouse_id,'_', '\_');
    
    % Get Unit Depth
    unit_depth = fetch1(data.Series(thisUnit), 'series_depth')-50-(25*(iunit-1));
    
    
    
    %%%%%%%%%%%
    %  CHIRP  %
    %%%%%%%%%%%
    
    if(flag.chirp == 1)
        try
            % Get key
            thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_chirp;
            sdfs = plotChirp(thisUnit, true, pathway);
            
            % Save the sdfs in a variable
            chirp_sdfs{iunit} = sdfs;
            
        catch
            display('skipped - chirp exception')
            
            % Save error message
            fprintf(stim_exceptions, '%s\t', 'Chirp');
            fprintf(stim_exceptions, '%s\t', 'Something went wrong');
            fprintf(stim_exceptions, '%s\t', mouse_id);
            fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
            fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
            fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            
            % Check for override
            if(plot_iunit > 0)
                break
            else
                continue;
            end
        end
        
        
        % Get monitor angle and elevation; used as standard
        [chirp_angle, chirp_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
        
        % Get Title
        supertitle(sprintf('Exp: %d, Chirp', thisUnit.exp_num));
        
        % Copy to layout
        figHandles = findall(gcf, 'Type', 'axes');
        newH = copyobj(figHandles(1), fh);   % stim
        set(newH, axPars, 'Position', [2.5 h-7.7 10 0.5], 'YTick', [], 'YColor','w', 'ticklength', [0.025, 0]);   % left, top, width, height
        
        newH = copyobj(figHandles(2), fh);   % text
        set(newH, axPars, 'Position', [4 h-2.2 10 1]);
        
        newH = copyobj(figHandles(3), fh);   % spikes
        set(newH, axPars, 'Position', [2.5 h-6.85 10 1.5], 'XTick', [], 'XColor','w','XTickLabel','');
        
        newH = copyobj(figHandles(4), fh);   % SDF
        set(newH, axPars, 'Position', [2.5 h-5 10 2]);
        
        close(gcf)
        
        % Add Colormap
        figure;
        colormap('jet')
        imagesc(mean(sdfs,1))
        set(gca, 'XTick', [],'YTick', []);
        
        % Copy to layout
        figHandles = findall(gcf, 'Type', 'axes');
        newH = copyobj(figHandles, fh);
        set(newH, axPars, 'Position', [2.5 h-7.1 10 0.2], 'XTick', [],'YTick', [], 'XColor','w','YColor','w','XTickLabel','');
        
        close(gcf)
        display '... Chirp done'
    end
    
    
    
    %%%%%%%%%
    % INFO  %
    %%%%%%%%%
    
    % Format ranksum
    ranksumf = sprintf('%.4f',thisUnit.ranksum);
    if(str2double(ranksumf) == 0)
        ranksumf = sprintf('1e%d',floor(log10(thisUnit.ranksum)));
    end
    
    % Add supertitle
    supertitle(sprintf('%s, SERIES: %d, UNIT: %d, DEPTH: %dum, P: %s, Qi: %.4f',latex_mouse_id, thisUnit.series_num, thisUnit.unit_id, unit_depth, ranksumf, thisUnit.qi));
    
    % Update key
    [thisUnit.clu_file_num, thisUnit.cluster_num] = fetchn(data.ClusterInfo(thisUnit), 'clu_file_num', 'cluster_num');
    
    % Get pathway to the clue files
    [~, ~, sortingfolder] = getPathTo(rmfield(thisUnit, {'clu_file_num', 'cluster_num'}), 'data');
    [~, basename, ~]=fileparts(sortingfolder);
    
    % Get the timestamps (for acorr)
    if(flag.autoc)
        clufile = [sortingfolder filesep basename '.clu.' num2str(thisUnit.clu_file_num)];
        clufid = fopen(clufile);
        cluids = fscanf(clufid,'%d');
        fclose(clufid);
        cluids = cluids(2:end);
        resfile = [sortingfolder filesep basename '.res.' num2str(thisUnit.clu_file_num)];
        resfid = fopen(resfile);
        clutimes = fscanf(resfid,'%d');
        fclose(resfid);
        timestamps = clutimes(cluids == thisUnit.cluster_num);
    end
    
    % Get avg spike waveforms
    if(flag.asw)
        spk = spkWaves(sortingfolder, thisUnit.clu_file_num, thisUnit.cluster_num);
        thisUnit.avgWave = mean(spk.waves,3);
        thisUnit.semWave = std(single(spk.waves),1,3)./sqrt(size(spk.waves,3));
        [thisUnit.aCorr, thisUnit.t] = PointCorrel(timestamps, timestamps, 1/1000*30000,30, 1, 30000, 'hz');
    end
    
    
    %%%%%%%%%%%
    %  AcTun  %
    %%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_cTun) && flag.ac)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_cTun);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_cTun(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING cTun: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_cTun(istim);
                plot(data.ContrastTuning(thisUnit));
                
                % Adjust graph
                title(sprintf('Exp: %d, Contrast Tuning', thisUnit.exp_num),'Interpreter', 'none');
                nBGHandle = get(gcf,'Children');
                set(nBGHandle(1),'Color', [1 1 1 ]);
                
                % Plot the best suited, depending on the angle and elevation,
                % and save the others
                if(istim == stim_num)
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
                    
                    % Copy to layout
                    figHandles = findall(gcf, 'Type', 'axes');
                    newH = copyobj(figHandles, fh);
                    set(newH, axPars, 'Position', [3 h-14.2 5.5 4]);
                    
                    close(gcf)
                    display '... cTun done'
                else
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num, istim)));
                    
                    close(gcf)
                    display '... Repeated cTun presentation';
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... AcTun skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'AcTun');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %  SparseNoise4Blank   %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_sparseNoise4Blank) && flag.sparse)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_sparseNoise4Blank);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_sparseNoise4Blank(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING SparseNoise4Blank: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_sparseNoise4Blank(istim);
                plot(data.SparseNoiseTuning(thisUnit));
                
                % Sort the figures
                for igraph = 1:4
                    
                    % close all fits
                    sGraphName = get(gcf,'Name');
                    if(strcmp(char(regexp(sGraphName,'fits','match','ignorecase')), 'fits') == 1)
                        close(gcf)
                        continue
                    end
                    
                    % Get the RF type & adjust graphs
                    if(strcmp(char(regexp(sGraphName,'-on-','match','ignorecase')), '-on-') == 1)
                        title(sprintf('Exp: %d, ON', thisUnit.exp_num));
                        sRFType = 'on';
                    else
                        title('OFF');
                        sRFType = 'off';
                    end
                    
                    % Plot the best suited, depending on the angle and
                    % elevation, and save the others
                    if(istim == stim_num)
                        % Save
                        print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%s.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num,sRFType)));
                        
                        % Copy to layout
                        figHandles = findall(gcf, 'Type', 'axes');
                        newH = copyobj(figHandles, fh);
                        
                        % Copy to layout accordingly to RF type
                        if(strcmp(sRFType,'on'))
                            set(newH, axPars, 'Position', [14 h-6.9 2.5 2.5]);    % ON Field
                        else
                            set(newH, axPars, 'Position', [17, h-6.9 2.5 2.5]);   % OFF Field
                        end
                        
                        display '... SparseNoise4Blank done'
                        
                        % Close the current figure
                        close(gcf)
                        
                    else
                        % Save
                        print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%s_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num,sRFType,istim)));
                        
                        display '... Repeated sparseNoise4Blank presentation'
                        
                        % Close the current figure
                        close(gcf)
                    end
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... SparseNoise4Blank skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'SparseNoise4Blank');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%
    %   AoriTun   %
    %%%%%%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_oriTun) && flag.aori)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_oriTun);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_oriTun(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING AoriTun: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_oriTun(istim);
                plot(data.OriTuning(thisUnit));
                xlabel('Direction (deg)');
                
                % Adjust graph
                title(sprintf('Exp: %d, Ori Tuning', thisUnit.exp_num),'Interpreter', 'none');
                nBGHandle = get(gcf,'Children');
                set(nBGHandle(1),'Color', [1 1 1 ]);
                
                % Plot the best suited, depending on the angle and elevation,
                % and save the others
                if(istim == stim_num)
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
                    
                    % Copy to layout
                    figHandles = findall(gcf, 'Type', 'axes');
                    newH = copyobj(figHandles, fh);
                    set(newH, axPars, 'Position', [3 h-20.7 5.5 4]);
                    
                    close(gcf)
                    display '... AoriTun done'
                else
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num, istim)));
                    
                    close(gcf)
                    display '... Repeated AoriTun presentation';
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... AoriTun skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'AoriTun');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    
    %%%%%%%%%%%%%%
    %   AtfTun   %
    %%%%%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_tfTun) && flag.atf)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_tfTun);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_tfTun(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING AtfTun: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_tfTun(istim);
                plot(data.TempFreqTuning(thisUnit));
                
                % Adjust graph
                title(sprintf('Exp: %d, Temp Freq Tuning', thisUnit.exp_num),'Interpreter', 'none');
                nBGHandle = get(gcf,'Children');
                set(nBGHandle(1),'Color', [1 1 1 ]);
                
                % Plot the best suited, depending on the angle and elevation,
                % and save the others
                if(istim == stim_num)
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
                    
                    % Copy to layout
                    figHandles = findall(gcf, 'Type', 'axes');
                    newH = copyobj(figHandles, fh);
                    set(newH, axPars, 'Position', [12.5 h-14.2 5.5 4]);
                    
                    close(gcf)
                    display '... AtfTun done'
                else
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num, istim)));
                    
                    close(gcf)
                    display '... Repeated AtfTun presentation';
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... AoriTun skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'AoriTun');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    %%%%%%%%%%%%%%%
    %   SizeTun   %
    %%%%%%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_sizeTun) && flag.size)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_sizeTun);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_sizeTun(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING sizeTun: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_sizeTun(istim);
                plot(data.SizeTuning(thisUnit));
                
                % Adjust graph
                title(sprintf('Exp: %d, Size Tuning', thisUnit.exp_num),'Interpreter', 'none');
                nBGHandle = get(gcf,'Children');
                set(nBGHandle(1),'Color', [1 1 1 ]);
                
                % Plot the best suited, depending on the angle and elevation,
                % and save the others
                if(istim == stim_num)
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
                    
                    % Copy to layout
                    figHandles = findall(gcf, 'Type', 'axes');
                    newH = copyobj(figHandles, fh);
                    set(newH, axPars, 'Position', [12.5 h-20.7 5.5 4]);
                    
                    close(gcf)
                    display '... sizeTun done'
                else
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num, istim)));
                    
                    close(gcf)
                    display '... Repeated sizeTun presentation';
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... SizeTun skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'SizeTun');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    
    
    
    %%%%%%%%%%%%%
    %   DsTun   %
    %%%%%%%%%%%%%
    
    if(~isempty(units_for_chirp_sorted(iunit).exp_num_DS) && flag.ds)
        
        % Has the same stimulus been repeated?
        [nStimNum,~] = size(units_for_chirp_sorted(iunit).exp_num_DS);
        
        % 0 = elevation & angle different from chirp, otherwise contains
        % the appropriate stimulus number
        stim_num = 0;
        
        % Choose stimulus with the same monitor angle & elevation as the
        % chirp stimulus
        if(nStimNum > 1)
            for istim = 1:nStimNum
                
                % Get angle and elevation for this stimulus
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_DS(istim);
                [this_angle, this_elevation] = fetchn(data.Experiments(thisUnit), 'exp_monitorangle', 'exp_monitorelevation');
                
                if(this_angle==chirp_angle && this_elevation==chirp_elevation)
                    stim_num = istim;
                    break
                end
            end
        else
            stim_num = 1;
        end
        
        % Parameters differ from the chirp stimulus, take the first one
        if(stim_num == 0)
            stim_num = 1;
            display '... WARNING DsTun: Different elevation or monitor angle';
        end
        
        for istim = 1:nStimNum
            try
                % Get key
                thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_DS(istim);
                dsTun(thisUnit);
                
                % Adjust graph
                title(sprintf('Exp: %d, DS Tuning', thisUnit.exp_num),'Interpreter', 'none');
                nBGHandle = get(gcf,'Children');
                set(nBGHandle(1),'Color', [1 1 1 ]);
                
                % Plot the best suited, depending on the angle and elevation,
                % and save the others
                if(istim == stim_num)
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
                    
                    % Copy to layout
                    figHandles = findall(gcf, 'Type', 'axes');
                    newH = copyobj(figHandles, fh);
                    set(newH, axPars, 'Position', [3 h-27 4 4]);
                    
                    close(gcf)
                    display '... DS done'
                else
                    % Save
                    print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d_%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num, istim)));
                    
                    close(gcf)
                    display '... Repeated DS presentation';
                end
            catch
                % Exception handling - Save exception error into as txt
                % file
                display('... DS skipped')
                
                % Save error message
                fprintf(stim_exceptions, '%s\t', 'DS');
                fprintf(stim_exceptions, '%s\t', 'Plot Error');
                fprintf(stim_exceptions, '%s\t', mouse_id);
                fprintf(stim_exceptions, '%d\t', thisUnit.series_num);
                fprintf(stim_exceptions, '%d\t', thisUnit.exp_num);
                fprintf(stim_exceptions, '%d\t\n', thisUnit.unit_id);
            end
        end
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Average Spike Wave %
    %%%%%%%%%%%%%%%%%%%%%%
    
    if(flag.asw == 1)
        
        % Plot
        figure
        plot(thisUnit.avgWave');
        title('Average Spike Wave');
        xlabel('Time (ms)');
        ylabel('\muV');
        
        % Save
        print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_ASW.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id)));
        
        % Copy to layout
        figHandles = findall(gcf, 'Type', 'axes');
        newH = copyobj(figHandles, fh);
        set(newH, axPars, 'Position', [2 h-32.7 6.5 4]);
        
        % Scale x-axis
        axisHandles = findall(fh, 'Type', 'axes');
        xt=get(axisHandles(1),'xticklabel');
        set(axisHandles(1),'xticklabel',num2str(round((str2double(xt)/30)*100)/100));
        
        close(gcf)
        display '... Average Spike Wave done'
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%
    % Autocorrelogram %
    %%%%%%%%%%%%%%%%%%%
    if(flag.autoc == 1)
        
        % Plot
        figure
        bar(thisUnit.t, thisUnit.aCorr, 'k');
        set(gca, axPars, 'xlim', [-20 20], 'ylim',[0 max(thisUnit.aCorr(thisUnit.t>=-20 & thisUnit.t <=20))]);
        title('Autocorrelogram');
        xlabel('Time lag (ms)');
        ylabel('Spikes/s');
        
        % Save
        print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_ACG.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id)));
        
        % Copy to layout
        figHandles = findall(gcf, 'Type', 'axes');
        newH = copyobj(figHandles, fh);
        set(newH, axPars, 'Position', [12 h-32.7 6.5 4]);
        
        close(gcf)
        display '... Autocorrelogram done'
    end
    
    
    % Save the layout
    if(flag.save == 1)
        set(gcf,'renderer','opengl');
        path = fullfile(pathway, 'Layout');
        if(exist(path, 'dir') ~= 7)
            mkdir(path)
        end
        
        export_fig(fullfile(path,sprintf('L_%s_s%02d_u%d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id)),'-nocrop','-eps');
        close(gcf)
        display(' ')
    end
    
    % Plot all or one unit; 0 = disable
    if(plot_iunit > 0 && numel(units_for_chirp_sorted) > 1)
        break
    end
end

% Save the mat files
if(flag.save == 1)
    path = fullfile(pathway, 'Data');
    if(exist(path, 'dir') ~= 7)
        mkdir(path)
    end
    save(fullfile(path, 'units_for_chirp_sorted.mat'),'units_for_chirp_sorted');
    save(fullfile(path, 'chirp_sdfs.mat'),'chirp_sdfs');
end










% AsfTun
%     if(~isempty(units_for_chirp_sorted(iunit).exp_num_AsfTun))
%         thisUnit.exp_num = units_for_chirp_sorted(iunit).exp_num_AsfTun;
%         plot(data.OriTuning(thisUnit));
%         print(gcf, '-depsc', fullfile(pathway, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, thisUnit.series_num, thisUnit.unit_id, thisUnit.exp_num)));
%
%         figHandles = findall(gcf, 'Type', 'axes');
%         newH = copyobj(figHandles, fh);
%         set(newH, axPars, 'Position', [12.5 18 5.5 4]);
%
%         close(gcf)
%         display '... AsfTun done'
%     end
%















