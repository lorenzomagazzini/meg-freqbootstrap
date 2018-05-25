
clear

addpath('/cubric/data/c1356674/ft_freqbootstrap')
addpath('/cubric/data/c1356674/ft_freqbootstrap/subfun')
addpath('/cubric/data/c1356674/ft_freqbootstrap/plotting')
data_path = '/cubric/scratch/c1356674/freqbootstrap_testing_20170301/data';
save_path = '/cubric/scratch/c1356674/freqbootstrap_testing_20170301';
dataset = fullfile(data_path, '230713-50_Standard_20170201_02.ds');


%% epoch the data

trial_timewin = [-1.5 1.5];

cfg = [];
cfg.dataset = dataset;
cfg.trialdef.eventtype  = 'Stim';
cfg.trialdef.prestim    = trial_timewin(1) * -1;
cfg.trialdef.poststim   = trial_timewin(2);
cfg_data = ft_definetrial(cfg);

% cfg_data.padding        = 2 * (cfg.trialdef.prestim + cfg.trialdef.poststim);
% cfg_data.padtype        = 'mirror';
cfg_data.continuous     = 'yes';
cfg_data.channel        = {'MEG'};
cfg_data.demean         = 'yes';
% cfg_data.hpfilter       = 'yes';
% cfg_data.hpfreq         = 1;
cfg_data.precision      = 'single';
data_raw = ft_preprocessing(cfg_data);

cd(save_path)
cd matfiles
save('data_raw.mat', '-v7.3', 'data_raw')

%% FFT of planar gradiometers

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.planarmethod = 'sincos';
cfg.neighbours   = ft_prepare_neighbours(cfg, data_raw);
data_planar = ft_megplanar(cfg, data_raw);

cfg = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.keeptrials  = 'yes';
cfg.pad         = 'nextpow2';
% cfg.polyremoval = 1;
% cfg.foilim      = [0 150];
% cfg.tapsmofrq   = [];
cfg.taper       = 'hanning';
data_planar_fft = ft_freqanalysis(cfg, data_planar);

cfg               = [];
cfg.combinemethod = 'sum';
data_fft = ft_combineplanar(cfg, data_planar_fft);

%% plot alpha power

% close all
% 
% cfg = [];
% cfg.parameter   = 'powspctrm';
% cfg.xlim        = [0 120];
% cfg.layout      = 'CTF275.lay';
% figure
% ft_multiplotER(cfg, data_fft)
% 
% cfg = [];
% cfg.parameter   = 'powspctrm';
% cfg.xlim        = [8 12];
% cfg.layout      = 'CTF275.lay';
% figure
% ft_topoplotER(cfg, data_fft)

%% select peak alpha channel

cfg = [];
cfg.keeptrials = 'no';
data_fft_avg = ft_freqdescriptives(cfg, data_fft);

foilim = [8 12];
freqs = data_fft_avg.freq;
numchan = length(data_fft_avg.label);

peakampl_arr = nan(numchan,1);
peakfreq_arr = nan(numchan,1);
for ichan = 1:numchan
    [peakampl_arr(ichan,1), peakfreq_arr(ichan,1)] = qcFindPeaksWithinFreqLims(data_fft_avg.powspctrm(ichan,:), freqs, foilim);
end

peakchannel_index = find(peakampl_arr == max(peakampl_arr));
peakchannel_label = data_fft_avg.label{peakchannel_index};

cfg = [];
cfg.channel = peakchannel_index;
data = ft_selectdata(cfg, data_fft)

cd(save_path)
cd matfiles
save('data_fft.mat', '-v7.3', '-struct', 'data')

%% bootstrapping

str = '_peakchan';

cfg = [];
cfg.parameter = 'powspctrm';
% cfg.operation = [];
cfg.foilim = [8 12];% [15 30];% 
cfg.findpeaks = 'yes';
cfg.findtroughs = 'no';
cfg.numboot = 10000;
cfg.winwidth = 1;
cfg.prctiter = 100;

outboot = ft_freqbootstrap(cfg, data);

cd(save_path)
cd matfiles
savename = ['outboot_' num2str(cfg.foilim(1)) '-' num2str(cfg.foilim(2)) 'Hz' str '.mat'];
save(savename, '-v7.3', '-struct', 'outboot')

%% plot spectra

xlim_spectra = [0 50];% [15 50];% 

close all

figure, set(gcf,'Color',[1 1 1])

hs1 = subplot(2,2,1);
hold on
kPlotShadedError(outboot.freq, outboot.spectra_bootmean, outboot.spectra_bootci95/2, [0.6 0.6 0.8])
plot(outboot.freq, outboot.spectra_bootmean, 'Color',[0 0 1], 'LineWidth',1)
plot(outboot.freq, outboot.spectra_bootmean-outboot.spectra_bootstdv, 'Color',[0 0 0.8], 'LineStyle','--')
plot(outboot.freq, outboot.spectra_bootmean+outboot.spectra_bootstdv, 'Color',[0 0 0.8], 'LineStyle','--')
xlim(xlim_spectra)
hl1 = plot([1 1]*outboot.maxfreq, get(gca,'YLim'), 'Color',[0 0 0.8], 'LineStyle','-.');
ylabel('Power (T^{2}/Hz)')
hleg1 = legend(hl1,char({'Peak Freq' [num2str(outboot.maxfreq,'%.2f') ' Hz']}));
hleg1.Box = 'off';

hs2 = subplot(2,2,2);
hold on
barh(outboot.maxampl_bootbins, outboot.maxampl_boothist, 'FaceColor',[0.6 0.6 0.8], 'EdgeColor',[0.6 0.6 0.8])
% xlim([0 outboot.cfg.numboot])
ylim(get(hs1,'YLim'))
xlabel('Summed iterations')

subplot(2,2,3)
hold on
plot(outboot.maxfreq_bootbins, outboot.maxfreq_boothist, 'Color',[0 0 1], 'LineWidth',1)
xlim(xlim_spectra)
ylim([0 outboot.cfg.numboot])
hl3 = plot([1 1]*outboot.maxfreq_bootmode, [0 outboot.cfg.numboot], 'Color',[0 0 0.8], 'LineStyle','-.');
xlabel('Frequency (Hz)')
ylabel('Summed iterations')
hleg3 = legend(hl3,char({'Boot Mode' [num2str(outboot.maxfreq_bootmode,'%.2f') ' Hz']}));
hleg3.Box = 'off';

cd(save_path)
cd figures
saveas(gcf, 'ft_freqbootstrap_example_alphapeakchan.png')
% saveas(gcf, 'ft_freqbootstrap_example_alphapeakchan_betafreq.png')

