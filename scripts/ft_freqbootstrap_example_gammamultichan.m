
clear

addpath('/cubric/data/c1356674/ft_freqbootstrap')
addpath('/cubric/data/c1356674/ft_freqbootstrap/subfun')
addpath('/cubric/data/c1356674/ft_freqbootstrap/plotting')
data_path = '/cubric/scratch/c1356674/freqbootstrap_testing_20170301/matfiles';


%% epoch the data

% trial_timewin = [-1.5 1.5];
% 
% cfg = [];
% cfg.dataset = dataset;
% cfg.trialdef.eventtype  = 'Stim';
% cfg.trialdef.prestim    = trial_timewin(1) * -1;
% cfg.trialdef.poststim   = trial_timewin(2);
% cfg_data = ft_definetrial(cfg);
% 
% % cfg_data.padding        = 2 * (cfg.trialdef.prestim + cfg.trialdef.poststim);
% % cfg_data.padtype        = 'mirror';
% cfg_data.continuous     = 'yes';
% cfg_data.channel        = {'MEG'};
% cfg_data.demean         = 'yes';
% % cfg_data.hpfilter       = 'yes';
% % cfg_data.hpfreq         = 1;
% cfg_data.precision      = 'single';
% data_raw = ft_preprocessing(cfg_data);

%% FFT of planar gradiometers

% cfg              = [];
% cfg.feedback     = 'no';
% cfg.method       = 'template';
% cfg.planarmethod = 'sincos';
% cfg.neighbours   = ft_prepare_neighbours(cfg, data_raw);
% data_planar = ft_megplanar(cfg, data_raw);
% 
% %separate conditions
% base_timewin = [-1.2 0];
% stim_timewin = [0.3 1.5];
% cfg = [];
% cfg.latency = base_timewin;
% data_planar_base = ft_selectdata(cfg, data_planar);
% cfg.latency = stim_timewin;
% data_planar_stim = ft_selectdata(cfg, data_planar);
% 
% cfg = [];
% cfg.method      = 'mtmfft';
% cfg.output      = 'pow';
% cfg.keeptrials  = 'yes';
% cfg.pad         = 'nextpow2';
% % cfg.polyremoval = 1;
% % cfg.foilim      = [0 150];
% % cfg.tapsmofrq   = [];
% cfg.taper       = 'hanning';
% data_planar_base_fft = ft_freqanalysis(cfg, data_planar_base);
% data_planar_stim_fft = ft_freqanalysis(cfg, data_planar_stim);
% 
% cfg               = [];
% cfg.combinemethod = 'sum';
% data_base_fft = ft_combineplanar(cfg, data_planar_base_fft);
% data_stim_fft = ft_combineplanar(cfg, data_planar_stim_fft);

%% plot gamma percentage change

% data_base_fft_avg = ft_freqdescriptives([], data_base_fft);
% data_stim_fft_avg = ft_freqdescriptives([], data_stim_fft);
% 
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.operation = '((x1-x2)./x2)*100';
% data_diff_fft = ft_math(cfg, data_stim_fft_avg, data_base_fft_avg);
% 
% close all
% 
% cfg = [];
% cfg.parameter   = 'powspctrm';
% cfg.xlim        = [40 70];
% cfg.layout      = 'CTF275.lay';
% figure
% ft_multiplotER(cfg, data_diff_fft)
% 
% cfg = [];
% cfg.parameter   = 'powspctrm';
% cfg.xlim        = [40 70];
% cfg.zlim        = 'maxabs';
% cfg.layout      = 'CTF275.lay';
% cfg.colormap   	= colormap('jet');
% figure
% ft_topoplotER(cfg, data_diff_fft)

%% bootstrapping

% str = '_multichandiff';
% 
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.operation = []; %use the default '((x1-x2)./x2)*100'
% cfg.foilim = [40 70];
% cfg.findpeaks = 'yes';
% cfg.findtroughs = 'no';
% cfg.numboot = 5000; %down to 5000 for computational reasons
% cfg.winwidth = 2;
% cfg.prctiter = 50; %50 percent iterations because gamma is noisier
% 
% outboot = ft_freqbootstrap(cfg, data_stim_fft, data_base_fft);
% 
% cd(data_path)
% savename = ['outboot_' num2str(cfg.foilim(1)) '-' num2str(cfg.foilim(2)) 'Hz' str '.mat'];
% save(savename, '-v7.3', '-struct', 'outboot')

%% plot topographies

str = '_multichandiff';
cfg.foilim = [40 70];
loadname = ['outboot_' num2str(cfg.foilim(1)) '-' num2str(cfg.foilim(2)) 'Hz' str '.mat'];
cd(data_path)
outboot = load(loadname);


close all
figure, set(gcf, 'Color',[1 1 1], 'Units','centimeters', 'Position',[5 17.5 40 7.5])

subplot(1,4,1)
cfg = [];
cfg.parameter = 'maxampl';
qcTopoplotSingleval(cfg, outboot)
title({'Peak amplitude' ''})

subplot(1,4,2)
cfg = [];
cfg.parameter = 'maxampl_bootmean';
qcTopoplotSingleval(cfg, outboot)
title({'Peak amplitude' '(bootstrap mean)'})

subplot(1,4,3)
cfg = [];
cfg.parameter = 'maxampl_bootmode';
qcTopoplotSingleval(cfg, outboot)
title({'Peak amplitude' '(bootstrap mode)'})

subplot(1,4,4)
cfg = [];
cfg.parameter = 'maxampl_bootstdv';
qcTopoplotSingleval(cfg, outboot)
title({'Peak amplitude' '(bootstrap SD)'})

cd(data_path)
cd ../figures
saveas(gcf, 'ft_freqbootstrap_example_gammamultichan_peakampl.png')


% close
figure, set(gcf, 'Color',[1 1 1], 'Units','centimeters', 'Position',[5 10 40 7.5])

subplot(1,4,1)
cfg = [];
cfg.parameter = 'maxfreq';
qcTopoplotSingleval(cfg, outboot)
title({'Peak frequency' ''})

subplot(1,4,2)
cfg = [];
cfg.parameter = 'maxfreq_bootmean';
qcTopoplotSingleval(cfg, outboot)
title({'Peak frequency' '(bootstrap mean)'})

subplot(1,4,3)
cfg = [];
cfg.parameter = 'maxfreq_bootmode';
qcTopoplotSingleval(cfg, outboot)
title({'Peak frequency' '(bootstrap mode)'})

subplot(1,4,4)
cfg = [];
cfg.parameter = 'maxfreq_prctiter';
% cfg.zlim = [0 100];
qcTopoplotSingleval(cfg, outboot)
title({'Percentage iterations' ['(' num2str(mode(outboot.cfg.winwidth),'%.1f') ' Hz around mode)']})

cd(data_path)
cd ../figures
saveas(gcf, 'ft_freqbootstrap_example_gammamultichan_peakfreq.png')




