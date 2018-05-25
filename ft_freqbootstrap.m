function [ outboot ] = ft_freqbootstrap( cfg, varargin )
%[ outboot ] = ft_freqbootstrap( cfg, ... )
%   
%   ft_freqbootstrap performs bootstrapping on spectral data over multiple 
%   trials. This function can be used to calculate measures of peak amplitude 
%   and peak frequency and to derive estimates of both central tendency and 
%   dispersion of these parameters across trials.
%   
%   Use as
%   	[outboot] = ft_freqbootstrap(cfg, freq)
%   or as
%     	[outboot] = ft_freqbootstrap(cfg, freq1, freq2)
%   
%   The input (freq) should be organised in a structure as obtained from the
%   FT_FREQANALYSIS function (with cfg.method='mtmfft', cfg.output='pow' and 
%   cfg.keeptrials='yes').
%   
%   The configuration should contain:
%     cfg.foilim      = [begin end], frequency band of interest.
%     cfg.findpeaks   = 'yes' or 'no', whether or not to find peaks in the
%                       power spectrum, within the frequency band specified 
%                       in cfg.foilim (default = 'yes').
%     cfg.findtroughs = 'yes' or 'no', whether or not to find troughs, i.e.
%                       negative peaks, within the frequency band specified 
%                       in cfg.foilim (default = 'no').
%     cfg.numboot     = scalar, number of bootstrap iterations (default = 10000).
%     cfg.winwidth    = scalar, width of the frequency window (in Hz) that
%                       will be used to calculate the percentage of bootstrap
%                       iterations around the bootstrap distribution mode
%                       (default = 2). Note that a value of 2 means +/- 1 Hz.                       
%     cfg.prctiter    = scalar, percentage of iterations that will be accommodated
%                       around the bootstrap distribution mode in order to
%                       calculate the necessary window width maxfreq_winwidth, 
%                       see below (default = 50).
%     cfg.parameter   = string indicating which field of the input data 
%                       structure will be used for the bootstrapping procedure 
%                       (default = 'powspctrm'). NOTE: 
%                       1) if a data structure with field powspctrm is passed on
%                       together with cfg.parameter = 'ampspctrm', then 
%                       the values in the powspctrm field will be square-rooted
%                       to derive spectral amplitude values. 
%                       2) if a data structure with field ampspctrm is passed on
%                       together with cfg.parameter = 'ampspctrm', then 
%                       the ampspctrm field will be assumed to contain spectral
%                       amplitude values and hence it will be left unchanged.
%   
%   If used as ... ft_freqbootstrap(cfg, freq1, freq2)
%     cfg.operation   = string indicating the mathematical operation that will 
%                       be performed on the (trial-averaged) field specified in 
%                       cfg.parameter (defalut = '((x1-x2)./x2)*100' ... i.e. 
%                       percentage change between freq1 and freq2). See also 
%                       help FT_MATH. The single-trial spectra will be averaged 
%                       (after each resampling iteration), before cfg.operation 
%                       is performed between x1 (for freq1) and x2 (for freq2). 
%   
%   The output of this function (outboot) is a stucture with fields:
%       maxampl_bootiter    = matrix (numchan x numboot), peak amplitude in each bootstrap iteration 
%       maxfreq_bootiter    = the same as maxampl_bootiter, but for peak frequency 
%       spectra_bootmean    = matrix (numchan x numfreq), spectrum averaged over bootstrap iterations 
%       spectra_bootstdv    = matrix (numchan x numfreq), standard deviation of the spectrum calculated across iterations 
%       spectra_bootci95    = matrix (numchan x numfreq), 95% confidence interval of the spectrum calculated across iterations 
%       maxampl             = vector (numchan x 1), peak amplitude in the original (non-bootstrapped) data 
%       maxampl_bootmean    = vector (numchan x 1), mean of the bootstrap peak amplitude distribution 
%       maxampl_bootmode    = vector (numchan x 1), mode of the bootstrap peak amplitude distribution 
%       maxampl_bootstdv    = vector (numchan x 1), SD of the bootstrap peak amplitude distribution
%       maxampl_boothist    = matrix (numchan x 100), number of iterations in each of the peak amplitude bins 
%       maxampl_bootbins    = matrix (numchan x 100), position of the peak amplitude bins 
%       maxfreq             = vector (numchan x 1), peak frequency in the original (non-bootstrapped) data 
%       maxfreq_bootmean    = the same as maxampl_bootmean, but for peak frequency 
%       maxfreq_bootmode    = the same as maxampl_bootmode, but for peak frequency 
%       maxfreq_bootstdv    = the same as maxampl_bootstdv, but for peak frequency 
%       maxfreq_boothist    = matrix (numchan x F, where F depends on cfg.foilim), the same as maxampl_boothist, but for peak frequency 
%       maxfreq_bootbins    = matrix (numchan x F, where F depends on cfg.foilim), the same as maxampl_bootbins, but for peak frequency 
%       maxfreq_prctiter    = vector (numchan x 1), percentage of iterations falling within the frequency window of width specified in cfg.winwidth, around maxfreq_bootmode 
%       maxfreq_winwidth    = vector (numchan x 1), width of the frequency window necessary to accommodate the percentage of iterations specified in cfg.prctiter, around maxfreq_bootmode
%   
%   If the configuration contains cfg.findtroughs = 'yes', then outboot
%   will contain the same fields as above, for the negative spectral peaks, 
%   named min* (e.g., minfreq) instead of max* (e.g., maxfreq).
%   
% 	To facilitate data-handling and distributed computing you can use 
%   cfg.inputfile     = ...
%   If you specify this, the input data will be read from a *.mat file on disk.
%   The mat files should contain only a single variable, corresponding to
%   the input structure. 
%   If the function is used with a single input data structure (freq),
%   cfg.inputfile should be a string.
%   If the function is used with two input structures (freq1 and freq2), 
%   cfg.inputfile should be a 2x1 cell array of strings, where the first
%   string corresponds to freq1 and the second string corresponds to freq2.


%% input cfg

%determine if the input is one or two spectra
switch nargin
    case 0
        eval('help ft_qc')
    case 1
        if ~isfield(cfg,'inputfile') || isempty(cfg.inputfile)
            error('this function requires data as input')
        else
            if ischar(cfg.inputfile)
                tmp = whos('-file',cfg.inputfile);
                tmp_name = tmp.name;
                fprintf('Loading file %s\n', cfg.inputfile)
                load(cfg.inputfile)
                singlespectra = eval(tmp_name);
                numcond = 1;
                clear tmp*
            elseif iscell(cfg.inputfile)
                tmp1 = whos('-file',cfg.inputfile{1});
                tmp1_name = tmp1.name;
                fprintf('Loading file %s\n', cfg.inputfile{1})
                load(cfg.inputfile{1})
                stimspectra = eval(tmp1_name);
                tmp2 = whos('-file',cfg.inputfile{2});
                tmp2_name = tmp2.name;
                fprintf('Loading file %s\n', cfg.inputfile{2})
                load(cfg.inputfile{2})
                basespectra = eval(tmp2_name);
                numcond = 2;
                clear tmp*
            end
        end
    case 2
        singlespectra = varargin{1};
        numcond = 1;
%         singlespectra = freq; %for debugging
    case 3
        stimspectra = varargin{1};
        basespectra = varargin{2};
        numcond = 2;
%         stimspectra = freq1; %for debugging
%         basespectra = freq2; %for debugging
    otherwise
        error('unexpected input')
end

%set defaults, check and throw errors
if ~isfield(cfg,'parameter') || isempty(cfg.parameter),                     cfg.parameter = 'powspctrm'; end
if ~isfield(cfg,'operation') || isempty(cfg.operation),                     cfg.operation = '((x1-x2)./x2)*100'; end
if numcond == 1,                                                           	cfg.operation = []; end
if ~isfield(cfg,'foilim') || isempty(cfg.foilim),                           error('input field cfg.foilim needs to be specified'); end
if ~isfield(cfg,'findpeaks') || isempty(cfg.findpeaks),                     cfg.findpeaks = 'yes'; end
if ~isfield(cfg,'findtroughs') || isempty(cfg.findtroughs),                	cfg.findtroughs = 'no'; end
if ~isfield(cfg,'numboot') || isempty(cfg.numboot),                         cfg.numboot = 10000; end
if ~isfield(cfg,'winwidth') || isempty(cfg.winwidth),                     	cfg.winwidth = 2; end
if ~isfield(cfg,'prctiter') || isempty(cfg.prctiter),                     	cfg.prctiter = 50; end

%compute spectral amplitude rather than power spectrum
if numcond == 1
    [singlespectra, cfg] = qcGetAmplitudeSpectra(cfg, singlespectra);
elseif numcond == 2
    [stimspectra, basespectra, cfg] = qcGetAmplitudeSpectra(cfg, stimspectra, basespectra);
end

%check cfg.parameter (ampspctrm, powspctrm)
if numcond == 1
    if ~isfield(singlespectra,cfg.parameter),                                       error(sprintf('The field ''%s'' is not present in the input data structure',cfg.parameter)); end
elseif numcond == 2
    if ~isfield(stimspectra,cfg.parameter) || ~isfield(basespectra,cfg.parameter),  error(sprintf('The field ''%s'' is not present in both input data structures',cfg.parameter)); end
end

%rename cfg input
parameter   = cfg.parameter;
operation   = cfg.operation;
foilim      = cfg.foilim;
findpeaks   = cfg.findpeaks;
findtrghs   = cfg.findtroughs;
numboot     = cfg.numboot;
winwidth    = cfg.winwidth;
prctiter    = cfg.prctiter;

%convert to logical
if strcmp(findpeaks,'yes'),     findpeaks = true;   else    findpeaks = false; end
if strcmp(findtrghs,'yes'),     findtrghs = true;   else    findtrghs = false; end

%get fields 'freq', 'label' and 'cfg.previous'
if numcond == 1
    
    label = singlespectra.label;
    freq = singlespectra.freq;
    
    if ~isfield(singlespectra,'cfg')
        previous = [];
    else
        previous = singlespectra.cfg;
    end
    
elseif numcond == 2
    
    label = stimspectra.label;
    if any(~strcmp(stimspectra.label, basespectra.label))
        error('this function requires freq1.label == freq2.label')
    end
    
    freq = stimspectra.freq;
    if any(stimspectra.freq ~= basespectra.freq)
        error('this function requires freq1.freq == freq2.freq')
    end
    
    if ~isfield(stimspectra,'cfg')
        previous{1} = [];
    else
        previous{1} = stimspectra.cfg;
    end
    if ~isfield(basespectra,'cfg')
        previous{2} = [];
    else
        previous{2} = basespectra.cfg;
    end
    
end


%% average

if numcond == 1
    %just average the single-trial spectra
    [avgspectrum, freqarray] = qcGetAverageSpectrum(singlespectra, parameter);
elseif numcond == 2
    %get spectrum of percentage change between stimulus and baseline
    [avgspectrum, freqarray] = qcGetAverageDiffSpectrum(stimspectra, basespectra, parameter, operation);
end

%% original data

%get number of channels
numchan = size(avgspectrum,1);
if numchan>1
%     warning('multi-channel functionality hasn''t been properly tested yet');
    fprintf('Bootstrapping will be repeated over %d channels\n', numchan)
end

%prepare the array to store the original peak parameters
maxpeakamplval = nan(numchan,1);
maxpeakfreqval = maxpeakamplval;
maxpeakfreqidx = maxpeakamplval;
mintrghamplval = maxpeakamplval;
mintrghfreqval = maxpeakamplval;
mintrghfreqidx = maxpeakamplval;

%loop over channels
for chan = 1:numchan
    powspctrm = avgspectrum(chan,:);
    %get peaks
    if findpeaks
        [maxpeakamplval(chan,1), maxpeakfreqval(chan,1), maxpeakfreqidx(chan,1)] = qcFindPeaksWithinFreqLims(powspctrm, freqarray, foilim);
    end
    %get troughs
    if findtrghs
        [mintrghamplval(chan,1), mintrghfreqval(chan,1), mintrghfreqidx(chan,1)] = qcFindTroughsWithinFreqLims(powspctrm, freqarray, foilim);
    end
end

%% bootstrap

if numcond == 1
    %number of trials
    numtrial = size(singlespectra.(parameter),1);
    numtrialstim = 0;
    numtrialbase = 0;
elseif numcond == 2
    %number of trials
    numtrial = 0;
    numtrialstim = size(stimspectra.(parameter),1);
    numtrialbase = size(basespectra.(parameter),1);
end

%get the resampled trial indeces for each bootstrap iteration
[bootindex, bootindexstim, bootindexbase] = qcGetResamplingIndex(numboot, numtrial, numtrialstim, numtrialbase);

%create cfg input to pass on to bootstrap function
cfg_boot = [];
cfg_boot.freq           = freq;
cfg_boot.label          = label;
cfg_boot.parameter      = parameter;
cfg_boot.operation      = operation;
cfg_boot.foilim         = foilim;
cfg_boot.findpeaks      = findpeaks;
cfg_boot.findtrghs      = findtrghs;
cfg_boot.numboot        = numboot;
cfg_boot.numchan        = numchan;
cfg_boot.numtrial       = numtrial;
cfg_boot.numtrialstim   = numtrialstim;
cfg_boot.numtrialbase   = numtrialbase;
cfg_boot.bootindex      = bootindex;
cfg_boot.bootindexstim  = bootindexstim;
cfg_boot.bootindexbase  = bootindexbase;

%do the bootstrapping
if numcond == 1
    [outboot, iterspectra] = qcBootstrapSingleSpectra(cfg_boot, singlespectra);
elseif numcond == 2
    [outboot, iterspectra] = qcBootstrapStimBaseSpectra(cfg_boot, stimspectra, basespectra);
end

%% get bootstrap summary measures

fprintf('Calculating bootstrap summary measures...\n')

%mean and standard deviation
if findpeaks
maxampl_bootmean = nanmean(outboot.maxampl_bootiter,2);
maxfreq_bootmean = nanmean(outboot.maxfreq_bootiter,2);
maxampl_bootstdv = nanstd(outboot.maxampl_bootiter,[],2);
maxfreq_bootstdv = nanstd(outboot.maxfreq_bootiter,[],2);
end
if findtrghs
minampl_bootmean = nanmean(outboot.minampl_bootiter,2);
minfreq_bootmean = nanmean(outboot.minfreq_bootiter,2);
minampl_bootstdv = nanstd(outboot.minampl_bootiter,[],2);
minfreq_bootstdv = nanstd(outboot.minfreq_bootiter,[],2);
end

%bootstrapped spectra mean and standard deviation
spectra_bootmean = mean(iterspectra,3);
spectra_bootstdv = nan(size(spectra_bootmean));
for chan = 1:numchan %doing this in a loop seems to prevent it from crashing...
    spectra_bootstdv(chan,:) = std(iterspectra(chan,:,:),[],3);
end

%define amplitude hist bin lims
foi = freqarray(nearest(freqarray, foilim(1), true, true) : nearest(freqarray, foilim(2), true, true));

%prepare arrays
numfreq = length(freqarray);
numbins = 100;
spectra_bootci95    = nan(numchan,numfreq);
amplcntrbins        = nan(numchan,numbins);
maxampl_bootbins    = amplcntrbins;
minampl_bootbins    = amplcntrbins;
maxampl_boothist    = amplcntrbins;
minampl_boothist    = amplcntrbins;
maxfreq_bootbins    = nan(numchan,length(foi));
minfreq_bootbins    = maxfreq_bootbins;
maxfreq_boothist    = maxfreq_bootbins;
minfreq_boothist    = maxfreq_bootbins;
maxampl_bootmode    = nan(numchan,1);
maxfreq_bootmode    = maxampl_bootmode;
minampl_bootmode    = maxampl_bootmode;
minfreq_bootmode    = maxampl_bootmode;

%loop over channels
for chan = 1:numchan
    
    tempspectra(1:numfreq,1:numboot) = squeeze(iterspectra(chan,:,:));
    
    %confidence interval of the bootstrapped spectra (95% C.I.)
    perc_vals = [2.5 97.5];
    perc = prctile(tempspectra, perc_vals, 2); %freq x perc
    spectra_bootci95(chan,1:numfreq) = abs(perc(:,2)-perc(:,1))'; %chan x freq
    
    %define amplitude hist bin lims (min and max ampl across iterations)
    perc_vals = [0 100];
    amplminmax = prctile(tempspectra, perc_vals, 2); %freq x perc
    amplcntrbins(chan,1:numbins) = linspace(min(amplminmax(:,1)), max(amplminmax(:,2)), numbins);
    
    clear tempspectra
    
    %derive mode from histogram distribution (don't use Matlab mode function..)
    if findpeaks
    [maxampl_boothist(chan,:), maxampl_bootbins(chan,:)] = hist(outboot.maxampl_bootiter(chan,:), amplcntrbins(chan,:));
    [histvalcnt, histvalidx] = max(maxampl_boothist(chan,:));
    maxampl_bootmode(chan,1) = maxampl_bootbins(chan,histvalidx);
%     figure, plot(maxampl_bootbins(chan,:), maxampl_boothist(chan,:), 'k'), hold on, plot(maxampl_bootmode(chan), histvalcnt, 'or'), plot([maxampl_bootmean(chan) maxampl_bootmean(chan)], [min(maxampl_boothist(chan,:)) max(maxampl_boothist(chan,:))], 'Color',[.5 .5 .5])
    clear histvalidx histvalcnt
    [maxfreq_boothist(chan,:), maxfreq_bootbins(chan,:)] = hist(outboot.maxfreq_bootiter(chan,:), foi);
    [histvalcnt, histvalidx] = max(maxfreq_boothist(chan,:));
    maxfreq_bootmode(chan,1) = maxfreq_bootbins(chan,histvalidx);
%     figure, plot(maxfreq_bootbins(chan,:), maxfreq_boothist(chan,:), 'k'), hold on, plot(maxfreq_bootmode(chan), histvalcnt, 'or'), plot([maxfreq_bootmean(chan) maxfreq_bootmean(chan)], [min(maxfreq_boothist(chan,:)) max(maxfreq_boothist(chan,:))], 'Color',[.5 .5 .5])
    clear histvalidx histvalcnt
    end
    if findtrghs
    [minampl_boothist(chan,:), minampl_bootbins(chan,:)] = hist(outboot.minampl_bootiter(chan,:), amplcntrbins(chan,:));
    [histvalcnt, histvalidx] = max(minampl_boothist(chan,:));
    minampl_bootmode(chan,1) = minampl_bootbins(chan,histvalidx);
%     figure, plot(minampl_bootbins(chan,:), minampl_boothist(chan,:), 'k'), hold on, plot(minampl_bootmode(chan), histvalcnt, 'or'), plot([minampl_bootmean(chan) minampl_bootmean(chan)], [min(minampl_boothist(chan,:)) max(minampl_boothist(chan,:))], 'Color',[.5 .5 .5])
    clear histvalidx histvalcnt
    [minfreq_boothist(chan,:), minfreq_bootbins(chan,:)] = hist(outboot.minfreq_bootiter(chan,:), foi);
    [histvalcnt, histvalidx] = max(minfreq_boothist(chan,:));
    minfreq_bootmode(chan,1) = minfreq_bootbins(chan,histvalidx);
%     figure, plot(minfreq_bootbins(chan,:), minfreq_boothist(chan,:), 'k'), hold on, plot(minfreq_bootmode(chan), histvalcnt, 'or'), plot([minfreq_bootmean(chan) minfreq_bootmean(chan)], [min(minfreq_boothist(chan,:)) max(minfreq_boothist(chan,:))], 'Color',[.5 .5 .5])
    clear histvalidx histvalcnt
    end
    
end

%% get the quality estimates

cfg_qc = [];
cfg_qc.winwidth    = winwidth; %requested width of frequency window around bootstrap mode
cfg_qc.prctiter    = prctiter; %requested percentage of iterations around the bootstrap mode
cfg_qc.numboot     = numboot; %requested number of bootstrap iterations

if findpeaks
    try
        maxfreq_bootiter = outboot.maxfreq_bootiter; %needed to calculate number of non-NaN iterations
        [maxfreq_prctiter, maxfreq_winwidth, cfg_qc] = qcGetPeakQualityEstimates(cfg_qc, maxfreq_boothist, maxfreq_bootbins, maxfreq_bootiter, maxfreq_bootmode);
    catch
        maxfreq_prctiter = NaN;
        maxfreq_winwidth = NaN;
        err = lasterror;
        outboot.maxfreq_errormsg = err.message;
    end
end
if findtrghs
    try
        minfreq_bootiter = outboot.minfreq_bootiter; %needed to calculate number of non-NaN iterations
        [minfreq_prctiter, minfreq_winwidth, cfg_qc] = qcGetPeakQualityEstimates(cfg_qc, minfreq_boothist, minfreq_bootbins, minfreq_bootiter, minfreq_bootmode);
    catch
        minfreq_prctiter = NaN;
        minfreq_winwidth = NaN;
        err = lasterror;
        outboot.minfreq_errormsg = err.message;
    end
end

%overwrite using actual winwidth used
cfg_boot.winwidth = cfg_qc.winwidth;
cfg_boot.prctiter = cfg_qc.prctiter;

%% append output fields

%bootstrapped spectra average, stdv and 95CI
outboot.spectra_bootmean = spectra_bootmean;
outboot.spectra_bootstdv = spectra_bootstdv;
outboot.spectra_bootci95 = spectra_bootci95;

%peaks (i.e. 'max')
if findpeaks

%max ampl orig % boot (mean, mode, stdv)
outboot.maxampl = maxpeakamplval;
outboot.maxampl_bootmean = maxampl_bootmean;
outboot.maxampl_bootmode = maxampl_bootmode;
outboot.maxampl_bootstdv = maxampl_bootstdv;

%max ampl hist & bins
outboot.maxampl_boothist = maxampl_boothist;
outboot.maxampl_bootbins = maxampl_bootbins;

%max freq orig % boot (mean, mode, stdv)
outboot.maxfreq = maxpeakfreqval;
outboot.maxfreq_bootmean = maxfreq_bootmean;
outboot.maxfreq_bootmode = maxfreq_bootmode;
outboot.maxfreq_bootstdv = maxfreq_bootstdv;

%max freq hist & bins
outboot.maxfreq_boothist = maxfreq_boothist;
outboot.maxfreq_bootbins = maxfreq_bootbins;

%QC
outboot.maxfreq_prctiter = maxfreq_prctiter;
outboot.maxfreq_winwidth = maxfreq_winwidth;

end

%troughs (i.e. 'min')
if findtrghs

%min ampl orig % boot (mean, mode, stdv)
outboot.minampl = mintrghamplval;
outboot.minampl_bootmean = minampl_bootmean;
outboot.minampl_bootmode = minampl_bootmode;
outboot.minampl_bootstdv = minampl_bootstdv;

%min ampl hist & bins
outboot.minampl_boothist = minampl_boothist;
outboot.minampl_bootbins = minampl_bootbins;

%min freq orig % boot (mean, mode, stdv)
outboot.minfreq = mintrghfreqval;
outboot.minfreq_bootmean = minfreq_bootmean;
outboot.minfreq_bootmode = minfreq_bootmode;
outboot.minfreq_bootstdv = minfreq_bootstdv;

%min freq hist & bins
outboot.minfreq_boothist = minfreq_boothist;
outboot.minfreq_bootbins = minfreq_bootbins;

%QC
outboot.minfreq_prctiter = minfreq_prctiter;
outboot.minfreq_winwidth = minfreq_winwidth;

end

%bookkeeping
outboot.cfg = cfg_boot;
outboot.cfg.previous = previous;


end

