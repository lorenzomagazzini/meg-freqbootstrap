function [ outboot, iterspectra ] = qcBootstrapSingleSpectra( cfg_boot, singlespectra )
%[ outboot, iterspectra ] = qcBootstrapSingleSpectra( cfg_boot, singlespectra )
%   Detailed explanation goes here

%get number of channels and number of iterations
numchan = cfg_boot.numchan;% size(singlespectra.powspctrm,2);
numboot = cfg_boot.numboot;

%prepare the array to store the bootstrapped peak parameters
if cfg_boot.findpeaks
maxpeakamplval = nan(numchan,numboot);
maxpeakfreqval = maxpeakamplval;
end

%prepare the array to store the bootstrapped trough parameters
if cfg_boot.findtrghs
mintrghamplval = nan(numchan,numboot);
mintrghfreqval = mintrghamplval;
end

%get frequency array
freqarray = singlespectra.freq;
numfreq = length(freqarray);

%prepare the array to store the single-iteration spectra
iterspectra = nan(numchan,numfreq,numboot);

%feedback
ft_progress('init', 'etf'); %comment out to maximise computational efficiency
fprintf('Bootstrapping...\n')
tic

%loop over iterations
for boot = 1:numboot
    ft_progress(boot/numboot, 'iteration %d of %d', boot, numboot); %comment out to maximise computational efficiency
    
    %prepare new chan x freq matrix
    bootavgspectrum = nan(numchan,numfreq);
    
    %loop over channels
    for chan = 1:numchan
        
        %average over resampled trials
        bootavgspectrum(chan,1:numfreq) = squeeze(mean(singlespectra.(cfg_boot.parameter)(cfg_boot.bootindex(boot,:),chan,:),1));
        
        %get peaks
        if cfg_boot.findpeaks
            [maxpeakamplval(chan,boot), maxpeakfreqval(chan,boot)] = qcFindPeaksWithinFreqLims(bootavgspectrum(chan,:), freqarray, cfg_boot.foilim);
        end
        
        %get troughs
        if cfg_boot.findtrghs
            [mintrghamplval(chan,boot), mintrghfreqval(chan,boot)] = qcFindTroughsWithinFreqLims(bootavgspectrum(chan,:), freqarray, cfg_boot.foilim);
        end
        
        %store single-iteration spectra
        iterspectra(chan,:,boot) = bootavgspectrum(chan,:);
        
    end %chan
    
end %boot
ft_progress('close') %comment out to maximise computational efficiency
fprintf('Finished.\n')
toc

%assign output
outboot         = struct;
outboot.freq	= cfg_boot.freq;
outboot.label	= cfg_boot.label;
if cfg_boot.findpeaks
    outboot.maxampl_bootiter = maxpeakamplval;
    outboot.maxfreq_bootiter = maxpeakfreqval;
end
if cfg_boot.findtrghs
    outboot.minampl_bootiter = mintrghamplval;
    outboot.minfreq_bootiter = mintrghfreqval;
end

end

