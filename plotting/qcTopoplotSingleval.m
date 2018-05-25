function [  ] = qcTopoplotSingleval( cfg, outboot )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

parameter = cfg.parameter;

colormapname = 'jet';

if isfield(cfg,'zlim') && ~isempty(cfg.zlim)
    zlim = cfg.zlim;
else
    zlim = [];
end

paramarray = {'maxfreq' 'maxfreq_bootmean' 'maxfreq_bootmode'};
if ismember(parameter,paramarray)
    zlim = outboot.cfg.foilim;% [floor(min(outboot.(parameter))) ceil(max(outboot.(parameter)))];
    colormapname = 'gray';
end

paramarray = {'maxfreq_prctiter'};
if ismember(parameter,paramarray)
    zlim = [0 100];
    colormapname = 'gray';
end

paramarray = {'maxampl' 'maxampl_bootmean' 'maxampl_bootmode' 'maxampl_bootstdv'};
if ismember(parameter,paramarray)
    if any(outboot.(parameter) > 0) && any(outboot.(parameter) < 0)
        zlim = 'maxabs';
    elseif any(outboot.(parameter) > 0)
        zlim = [0 max(outboot.(parameter))];
    elseif any(outboot.(parameter) < 0)
        zlim = [min(outboot.(parameter)) 0];
    end
end


ft_data = struct;
ft_data.dimord      = 'chan_freq';
ft_data.label       = outboot.label;
ft_data.freq        = 1;
ft_data.powspctrm   = outboot.(parameter);

cfg = [];
cfg.zlim        = zlim;
cfg.layout      = 'CTF275.lay';
cfg.comment     = 'no';
cfg.interpolatenan	= 'yes';
cfg.colorbar    = 'yes';
cfg.gridscale   = 200;
cfg.colormap    = colormap(colormapname);
% figure
ft_topoplotER(cfg, ft_data)


end

