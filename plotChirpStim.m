%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HARD-CODED CHIRP STIMULUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The function generates variables for plotting the chirp stimulus
%
% @t     ... timebase of the stimulus
% @y     ... y-values of the stimulus
% @cumTs ... timing of the different stimulus parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y, cumTs] = plotChirpStim()

% VARS
ChirpDuration = 8;                   % Time (s) of rising/falling chirp phase
ChirpMaxFreq  = 8;                   % Peak frequency of chirp (Hz)
IntensityFrequency = 2;               % freq at which intensity is modulated

% FREQ
Fduration     = 0.017;                % Single Frame duration (s) -  ADJUST DEPENDING ON MONITOR
Fduration_ms  = 17.0;                 % Single Frame duration (ms) - ADJUST DEPENDING ON MONITOR
KK            = ChirpMaxFreq / ChirpDuration; % acceleration in Hz / s
Npoints       = round(ChirpDuration / Fduration);
freq_stim = nan(1,Npoints);
for it = 1 : Npoints
    Time = it*Fduration; % in s
    freq_stim(it) = sin(3.141 * KK * Time*Time );
end

% CONTRAST
KK2           = IntensityFrequency;
c_stim = nan(1,Npoints);
for it = 1 : Npoints
    Time = it*Fduration;  % in s
    %RampIntensity = int(127*Time/ChirpDuration);
    RampIntensity = Time/ChirpDuration;
    c_stim(it) = sin(2*3.141 * KK2 * Time) * RampIntensity;      % gives flicker at 1/2 max Chirp freq
end


SteadyOFF     = 3.00;                 % Time (s) of Light OFF at beginning at end of stimulus
SteadyOFF2    = 2.00;
SteadyON      = 3.00;                 % Time (s) of Light 100% ON before and after chirp
SteadyMID     = 2.00;                 % Time (s) of Light at 50% for steps
SteadyMID2    = 1.00;                 % Time (s) of Light at 50% between flickers

% Define separate chirp parts
cumTs = cumsum([SteadyOFF2, SteadyON, SteadyOFF, SteadyMID, ChirpDuration, SteadyMID, ChirpDuration, SteadyMID, SteadyOFF2]);

t = (0:Fduration:cumTs(end)); % t in s
y = nan(size(t));
y(t>0 & t <=cumTs(1))        = -1;    % 2s_OFF    2
y(t>cumTs(1) & t <=cumTs(2)) = 1;     % 3s_ON     5
y(t>cumTs(2) & t <=cumTs(3)) = -1;    % 3s_OFF    8
y(t>cumTs(3) & t <=cumTs(4)) = 0;     % 1s_MID    9
y(t>cumTs(4) & t <=cumTs(5)) = freq_stim(1:nnz(t>cumTs(4) & t <=cumTs(5)));
y(t>cumTs(5) & t <=cumTs(6)) = 0;     % 2s_STEADYMID
y(t>cumTs(6) & t <=cumTs(7)) = c_stim(1:nnz(t>cumTs(6) & t <=cumTs(7)));
y(t>cumTs(7) & t <=cumTs(8)) = 0; % 2s_STEADYMID
y(t>cumTs(8) & t <=cumTs(9)) = -1; % 2s_STEADYOFF
end