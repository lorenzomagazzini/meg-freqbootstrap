function [ outboot, iterspectra ] = qcBootstrapStimBaseSpectra( cfg_boot, stimspectra, basespectra )
%[ outboot, iterspectra ] = qcBootstrapStimBaseSpectra( cfg_boot, stimspectra, basespectra )
%   Detailed explanation goes here

%TODO - check if basespectra.label == stimspectra.label (needed?)

%get number of channels
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
freqarray = stimspectra.freq;
numfreq = length(freqarray);

%prepare the array to store the single-iteration spectra
iterspectra = nan(numchan,numfreq,numboot);

%redefine operations string
operation = strrep(strrep(cfg_boot.operation,'x1','bootavgstimspectra(chan,1:numfreq)'),'x2','bootavgbasespectra(chan,1:numfreq)');

%feedback
ft_progress('init', 'etf'); %comment out to maximise computational efficiency
fprintf('Bootstrapping...\n')
tic

%loop over iterations
for boot = 1:numboot
    ft_progress(boot/numboot, 'iteration %d of %d', boot, numboot);
    
    %prepare new chan x freq matrices
    bootavgspectrum = nan(numchan,numfreq);
    bootavgstimspectra = bootavgspectrum;
    bootavgbasespectra = bootavgspectrum;
    
    %loop over channels
    for chan = 1:numchan
        
        %stimulus and baseline spectra resampling (avoid ft functions feedback)
        bootavgstimspectra(chan,1:numfreq) = squeeze(mean(stimspectra.(cfg_boot.parameter)(cfg_boot.bootindexstim(boot,:),chan,:),1));
        bootavgbasespectra(chan,1:numfreq) = squeeze(mean(basespectra.(cfg_boot.parameter)(cfg_boot.bootindexbase(boot,:),chan,:),1));
        
        %get percentage change spectrum
        bootavgspectrum(chan,1:numfreq) = eval(operation);
        
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

