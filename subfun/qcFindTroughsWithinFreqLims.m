function [ mintrghamplval, mintrghfreqval, mintrghfreqidx ] = qcFindTroughsWithinFreqLims( powspctrm, freqarray, foilim )
%[ mintrghamplval, mintrghfreqval, mintrghfreqidx ] = qcFindTroughsWithinFreqLims( powspctrm, freqarray, foilim )
%   
%   Input:
%   powspctrm:  1xF array of power values
%   freqarray:  1xF array of frequencies
%   foilim:     1x2 array indicating the lower and upper frequency limits
%   
%   Output:
%   mintrghamplval: amplitude of the trough with lowest amplitude within
%                   the specified frequency limits (foilim)
%   mintrghfreqval: frequency of the trough with lowest amplitude within
%                   the specified frequency limits (foilim)
%   mintrghfreqidx: index to the frequency corresponding to mintrghfreqval 
%                   in freqarray
%   

%define freq lims
f1_idx = nearest(freqarray, foilim(1), true, true);
% f1_val = freqarray(f1_idx);
f2_idx = nearest(freqarray, foilim(2), true, true);
% f2_val = freqarray(f2_idx);

%search for TROUGHS
%%%%% find troughs across the whole spectrum
alltrghs_idx = find(diff( diff(powspctrm)<0 ) <0) +1;
if isempty(alltrghs_idx),
    %%%%% if no peaks, then set as empty
    limtrghs_idx = [];
else
    %%%%% find troughs within the specified frequency range
    limtrghs_idx = alltrghs_idx(alltrghs_idx>=f1_idx & alltrghs_idx<=f2_idx);
end
if isempty(limtrghs_idx),
    %%%%% if no peaks, then set as NaN
    mintrghamplval = NaN;
    mintrghfreqval = NaN;
    mintrghfreqidx = NaN;
else
    %%%%% get MIN TROUGH amplitude and frequency
    [mintrghamplval, tmp_idx] = min(powspctrm(limtrghs_idx));
    mintrghfreqidx = limtrghs_idx(tmp_idx);
    mintrghfreqval = freqarray(mintrghfreqidx);
    clear tmp_idx
end

% %plot all troughs
% close all
% figure
% plot(powspctrm)
% hold on
% plot(alltrghs_idx, powspctrm(alltrghs_idx), 'ob')

% %plot troughs (within frequency range)
% close all
% figure
% plot(powspctrm)
% hold on
% plot(limtrghs_idx, powspctrm(limtrghs_idx), 'ob')

% %plot min trough
% close all
% figure
% plot(powspctrm)
% hold on
% plot(mintrghfreqidx, mintrghamplval, 'ob')

end

