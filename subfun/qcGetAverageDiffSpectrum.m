function [ avgspectrum, freqarray ] = qcGetAverageDiffSpectrum( stimspectra, basespectra, parameter, operation )
%[ avgspectrum, freqarray ] = qcGetAverageDiffSpectrum( stimspectra, basespectra, parameter, operation )
%   
%   Input:
%   stimspectra:    frequency structure, as obtained with ft_freqanalysis,
%                   containing the spectral data for the stimulus epoch
%   basespectra:    frequency structure, as obtained with ft_freqanalysis,
%                   containing the spectral data for the baseline epoch
%   parameter:      string, pointing to the structure field on which to
%                   apply the ft_math operation (...note that the field
%                   specified in parameter is renamed 'powspctrm', inside 
%                   the function)
%   operation:      ft_math operation between stimspectra and basespectra
%                   (x1 and x2, respectively)
%   
%   Output:
%   avgspectrum:        1xF array of power values, computed as a resulted
%                       of the ft_math 'operation' between stimspectra and 
%                       basespectra, applied to the field specified in
%                       'parameter'
%   freqarray:        	1xF array of frequencies for which the values in 
%                       avgspectrum were calculated
%                       

%TODO - check if basespectra.freq == stimspectra.freq (needed?)

%check dimord
if ~strcmp(stimspectra.dimord,'rpt_chan_freq')
    error('Wrong input data, stimspectra.dimord should be ''rpt_chan_freq''')
end
if ~strcmp(basespectra.dimord,'rpt_chan_freq')
    error('Wrong input data, basespectra.dimord should be ''rpt_chan_freq''')
end

%replace powspctrm field with field of interest
stimspectra.powspctrm = stimspectra.(parameter);
basespectra.powspctrm = basespectra.(parameter);

%average baseline and stimulus power spectra
cfg_avg = [];
cfg_avg.keeptrials = 'no';
stimspectra_avg = ft_freqdescriptives(cfg_avg, stimspectra);
basespectra_avg = ft_freqdescriptives(cfg_avg, basespectra);
clear cfg_avg

%average percentage change power spectrum
cfg_diff = [];
cfg_diff.parameter = 'powspctrm'; %(see field 'powspctrm' above)
cfg_diff.operation = operation;
diffspectra_avg = ft_math(cfg_diff, stimspectra_avg, basespectra_avg);
clear cfg_diff

%assign output
avgspectrum = diffspectra_avg.powspctrm;
freqarray = diffspectra_avg.freq;

end

