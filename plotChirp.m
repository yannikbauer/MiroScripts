%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT CHIRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots the chirp stimulus
%
% params:
% @key ... struct containing the pointer to the exp
% @printFlag ... 1=graph will be displayed
% @figPath ... path for saving the figure
%
% returns:
% @sdfs ... returns the spike density
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sdfs] = plotChirp(key, printFlag, figPath)

if nargin < 3
    figPath = pwd;
end

if nargin < 2
    printFlag = false;
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(data.Units(key),2,[], 10, [], 'extra');
[sdfs, ~] = plot(data.Units(key),2,[], 25, [], [], 'extra');

% Generates variables for the chirp stimulus
[t, y, cumTs] = plotChirpStim();




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT PARAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mouse_id = fetch1(data.Mice(key), 'mouse_id');

[figPars, axPars] = setPlotPars('slide');
axPars.TickLength = [0.025 0.025];
axPars.Units = 'centimeters';
myPos = [10 10 14 8.5];

figHandles = findall(gcf,'Type','figure');
axHandles = findall(gcf, 'Type', 'axes');
set(axHandles, axPars);

lineHandles = findall(axHandles, 'LineWidth', 1);
set(lineHandles, 'LineWidth', 0.5);

lineHandles = findall(axHandles, 'LineWidth', 2);
set(lineHandles, 'LineWidth', 0.5);

% remove the error
lineHandles = findall(axHandles, 'LineStyle', '--');
delete(lineHandles);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CHIRP & STIM FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the figure with both the chirp and the timeline
fh = figure(figPars, 'Position', myPos, 'PaperPosition', myPos);
clf;

spikeH = copyobj(axHandles(2),fh);
set(spikeH, axPars, 'Position', [1.5 2*(1+myPos(4)/4) myPos(3)-2 myPos(4)/4], ...
    'XTickLabel', [], 'Visible','off');

sdfH = copyobj(axHandles(1),fh);
set(sdfH, axPars, 'Position', [1.5 1+myPos(4)/4 myPos(3)-2 myPos(4)/4], ...
    'XTickLabel', [], 'XColor', 'w');
ylabel('Spikes/s');

ah = axes(axPars, 'Position', [1.5 1.5 myPos(3)-2 myPos(4)/12]);
plot(t,y, '-k');
set(ah, axPars, 'YLim', [-1.5 1.5], 'YColor', 'w')

set([ah sdfH spikeH], 'XLim', [0 cumTs(end)]);
xlabel('Time (s)');

figName = sprintf('%s_s%02d_e%02d_u%d', mouse_id, key.series_num, key.exp_num, key.unit_id);
set(gcf,'Name',figName);

close(figHandles);

% Save data
if printFlag
    print(gcf, '-deps', fullfile(figPath, sprintf('%s_s%02d_u%d_e%02d.eps', mouse_id, key.series_num, key.unit_id, key.exp_num)));
end
end



















