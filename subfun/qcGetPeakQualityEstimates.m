function [ prctiter_withinspecifiedwidth, winwidth_forspecifiedperciter, cfg ] = qcGetPeakQualityEstimates( cfg, freqhist, freqbins, freqiter, freqmode )
%[ prctiter_withinspecifiedwidth, winwidth_forspecifiedperciter, cfg ] = qcGetPeakQualityEstimates( cfg, freqhist, freqbins, freqiter, freqmode )
%   
%   Input:
%   cfg, structure with fields:
%      'prctiter':  requested percentage of iterations
%      'winwidth':  requested frequency window width
%      'numboot':   total number of number of bootstrap iterations
%   freqhist:       bootstrap (histogram) peak frequency distribution
%   freqbins:       frequencies of the bootstrap histogram distribution (in Hz)
%   freqmode:       mode of the bootstrap histogram distribution (in Hz)
%   
%   Output:
%   

% %debugging:
% cfg = cfg_qc;
% freqhist = maxfreq_boothist;
% freqbins = maxfreq_bootbins;
% freqiter = maxfreq_bootiter;
% freqmode = maxfreq_bootmode;
% cfg = cfg_qc;
% freqhist = minfreq_boothist;
% freqbins = minfreq_bootbins;
% freqiter = minfreq_bootiter;
% freqmode = minfreq_bootmode;


%warning
fprintf('Assuming equally-spaced frequency bins...\n') %see below

%requested percentage of iterations around bootstrap mode
prctiter = cfg.prctiter;

%number of channels
numchan = size(freqhist,1);

if any(size(cfg.winwidth) ~= [1 1]), error('cfg.winwidth should be a single scalar value'); end


%prepare output arrays
prctiter_withinspecifiedwidth = nan(numchan,1);
winwidth_forspecifiedperciter = nan(numchan,1);

%print feedback
tmp_deltafreq = unique(diff(freqbins(:,1:2)'));
tmp_freqarray = 0:tmp_deltafreq:tmp_deltafreq*1000;
tmp_halfwinwidth = tmp_freqarray(nearest(tmp_freqarray,cfg.winwidth/2));
fprintf('Calculating percentage of iterations within +/- %.2f Hz.\n', tmp_halfwinwidth)
clear tmp*

%repeat winwidth across channels
cfg.winwidth = repmat(cfg.winwidth,numchan,1);

for chan = 1:numchan
    
    %requested width of frequency window around bootstrap mode
    desired_winwidth = cfg.winwidth(chan,1);
    
    %maximum winwidth allowed
    maxhalfwinwidth = min([freqbins(chan,end)-freqmode(chan) freqmode(chan)-freqbins(chan,1)]); %distance between mode and closest frequency edge
    maxwinwidth = maxhalfwinwidth*2;
    
    if cfg.winwidth(chan,1) > ceil(maxwinwidth)
        
        warning(sprintf('The window width requested (%.2f Hz) is too large. The maximum width allowed for channel %d is %.2f Hz. Skipping...', cfg.winwidth(chan,1), chan, maxwinwidth));
        
        cfg.winwidth(chan,1) = NaN;
        prctiter_withinspecifiedwidth(chan,1) = NaN;
        
    else
        
        %frequency resolution (freq bins are assumed to be equally spaced)
        deltafreq = diff(freqbins(chan,1:2));
        
        %window half width (needs to be defined based on deltafreq steps)
        freqarray = 0:deltafreq:maxwinwidth;
        halfwinwidth = freqarray(nearest(freqarray,desired_winwidth/2));
        
        %actual window width to be used (overwrite based on delta freq)
        cfg.winwidth(chan,1) = halfwinwidth*2;
        
        
        %search window
        winsearchlim = [freqmode(chan)-halfwinwidth, freqmode(chan)+halfwinwidth];
        winsearchidx = nearest(freqbins(chan,:),winsearchlim(1)) : nearest(freqbins(chan,:),winsearchlim(2));
        
        %n. of NaN iterations and n. of non-NaN iterations
        numnaniter = sum(isnan(freqiter(chan,:))); 
        nonnaniter = sum(freqhist(chan,:));
        if nonnaniter ~= cfg.numboot-numnaniter, error('number of NaN iterations doesn''t match'); end
        
        
        %get percentage of iterations within specified window around bootstrap mode
        prctiter_withinspecifiedwidth(chan,1) = sum(freqhist(chan,winsearchidx)) / nonnaniter * 100;
        
        
    end %if...else
    
end %for...chan


%print feedback
fprintf('Calculating frequency window width necessary to accommodate %d%% of iterations.\n', prctiter)

%repeat prctiter across channels
cfg.prctiter = repmat(cfg.prctiter,numchan,1);

for chan = 1:numchan
    
    for deltaindx = 0:floor(length(freqbins)/2)
        
        clear tmp*
        tmp_winsearchlim = [freqmode(chan)-(deltaindx*deltafreq) freqmode(chan)+(deltaindx*deltafreq)];
        tmp_winsearchidx = nearest(freqbins(chan,:),tmp_winsearchlim(1)) : nearest(freqbins(chan,:),tmp_winsearchlim(2));
        
        %get percentage of iteration at current window width
        tmp_prctiter = sum(freqhist(chan,tmp_winsearchidx) / nonnaniter * 100);
        
        %check if requested percentage of iterations is exceeded, if so abort
        if tmp_prctiter >= cfg.prctiter(chan,1), break, end
        
    end
    
    
    %get window width accommodating specified percentage iterations around mode
    winwidth_forspecifiedperciter(chan,1) = freqbins(tmp_winsearchidx(end)) - freqbins(tmp_winsearchidx(1));
    
    %the following has been commented because an odd number of freq bins
    %may need to be used to go up to 100% iterations
    %     if winwidth_forspecifiedperciter(chan,1) ~= deltaindx*deltafreq*2, error('delta frequency doesn''t match'); end
    
    
end %for...chan

end

