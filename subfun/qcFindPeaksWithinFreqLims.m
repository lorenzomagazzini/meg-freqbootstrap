function [ maxpeakamplval, maxpeakfreqval, maxpeakfreqidx ] = qcFindPeaksWithinFreqLims( powspctrm, freqarray, foilim )
%[ maxpeakamplval, maxpeakfreqval, maxpeakfreqidx ] = qcFindPeaksWithinFreqLims( powspctrm, freqarray, foilim )
%   
%   Input:
%   powspctrm:  1xF array of power values
%   freqarray:	1xF array of frequencies
%   foilim:     1x2 array indicating the lower and upper frequency limits
%   
%   Output:
%   maxpeakamplval: amplitude of the peak with greatest amplitude within
%                   the specified frequency limits (foilim)
%   maxpeakfreqval: frequency of the peak with greatest amplitude within
%                   the specified frequency limits (foilim)
%   maxpeakfreqidx: index to the frequency corresponding to maxpeakfreqval 
%                   in freqarray
%   

%define freq lims
f1_idx = nearest(freqarray, foilim(1), true, true);
% f1_val = freqarray(f1_idx);
f2_idx = nearest(freqarray, foilim(2), true, true);
% f2_val = freqarray(f2_idx);

%search for PEAKS
%%%%% find peaks across the whole spectrum
allpeaks_idx = find(diff( diff(powspctrm)>0 ) <0) + 1;
if isempty(allpeaks_idx),
    %%%%% if no peaks, then set as empty
    limpeaks_idx = [];
else
    %%%%% find peaks within the specified frequency range
    limpeaks_idx = allpeaks_idx(allpeaks_idx>=f1_idx & allpeaks_idx<=f2_idx);
end
if isempty(limpeaks_idx),
    %%%%% if no peaks, then set as NaN
    maxpeakamplval = NaN;
    maxpeakfreqval = NaN;
    maxpeakfreqidx = NaN;
else
    %%%%% get MAX PEAK amplitude and frequency
    [maxpeakamplval, tmp_idx] = max(powspctrm(limpeaks_idx));
    maxpeakfreqidx = limpeaks_idx(tmp_idx);
    maxpeakfreqval = freqarray(maxpeakfreqidx);
    clear tmp_idx
end

% %plot all peaks
% close all
% figure
% plot(powspctrm)
% hold on
% plot(allpeaks_idx, powspctrm(allpeaks_idx), 'or')

% %plot peaks (within frequency range)
% close all
% figure
% plot(powspctrm)
% hold on
% plot(limpeaks_idx, powspctrm(limpeaks_idx), 'or')

% %plot max peak
% close all
% figure
% plot(powspctrm)
% hold on
% plot(maxpeakfreqidx, maxpeakamplval, 'or')

end

