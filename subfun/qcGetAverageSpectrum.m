function [ avgspectrum, freqarray ] = qcGetAverageSpectrum( singlespectra, parameter )
%[ avgspectrum, freqarray ] = qcGetAverageSpectrum( singlespectra, parameter )
%   
%   Input:
%   singlespectra:	frequency structure, as obtained with ft_freqanalysis,
%                   containing the single-trial spectra
%   parameter:      string, pointing to the field to average (...note that 
%                   the field specified in parameter renamed 'powspctrm', 
%                   inside the function)
%   
%   Output:
%   avgspectrum:    1xF array of power values, computed by averaging over  
%                   trials
%   freqarray:      1xF array of frequencies for which the values in 
%                   avgspectrum were calculated
%                       

%TODO - check if basespectra.freq == stimspectra.freq (needed?)

%check dimord
if ~strcmp(singlespectra.dimord,'rpt_chan_freq')
    error('Wrong input data, singlespectra.dimord should be ''rpt_chan_freq''')
end

%replace powspctrm field with field of interest
singlespectra.powspctrm = singlespectra.(parameter);

%average power spectra
cfg_avg = [];
cfg_avg.keeptrials = 'no';
cfg_avg.feedback = 'no';
spectra_avg = ft_freqdescriptives(cfg_avg, singlespectra);
clear cfg_avg

%assign output
avgspectrum = spectra_avg.powspctrm;
freqarray = spectra_avg.freq;

end

