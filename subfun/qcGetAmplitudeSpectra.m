function [ varargout ] = qcGetAmplitudeSpectra( cfg, varargin )
%[ varargout ] = qcGetAmplitudeSpectra( cfg, varargin )
%   Detailed explanation goes here

switch nargin
    case 2
        singlespectra = varargin{1};
        numcond = 1;
        fprintf('Input is one freq structure.\n')
    case 3
        stimspectra = varargin{1};
        basespectra = varargin{2};
        numcond = 2;
        fprintf('Input is two freq structures.\n')
    otherwise
        error('please check the use of this function')
end

if numcond == 1
    
    if strcmp(cfg.parameter,'ampspctrm') && isfield(singlespectra,'ampspctrm')
        fprintf('The field ''%s'' is assumed to consist of amplitude, rather than power values. Not doing anything to it...\n', cfg.parameter)
        
        if isfield(singlespectra,'powspctrm'), singlespectra = rmfield(singlespectra,'powspctrm'); fprintf('Removing powspctrm field...\n'); end
        
    elseif strcmp(cfg.parameter,'ampspctrm') && isfield(singlespectra,'powspctrm')
        fprintf('The field ''%s'' will be calculated by square-rooting the values in ''%s'', which are assumed to be power values.\n', cfg.parameter,'powspctrm')
        for trial = 1:size(singlespectra.powspctrm,1)
            for chan = 1:size(singlespectra.powspctrm,2)
                singlespectra.ampspctrm(trial,chan,:) = sqrt(squeeze(singlespectra.powspctrm(trial,chan,:)));
            end
        end
        
        singlespectra = rmfield(singlespectra,'powspctrm'); fprintf('Removing powspctrm field...\n');
        
    elseif strcmp(cfg.parameter,'ampspctrm')
        error('The cfg option parameter = ''%s'' can only be specified if either ''%s'' or ''%s'' is a field of the input data structure', cfg.parameter,'ampspctrm','powspctrm')
        
    end
    
    varargout{1} = singlespectra;% rmfield(singlespectra,'powspctrm');
    varargout{2} = cfg;
    
elseif numcond == 2
    
    if strcmp(cfg.parameter,'ampspctrm') && isfield(stimspectra,'ampspctrm') && isfield(basespectra,'ampspctrm')
        fprintf('The field ''%s'' is assumed to consist of amplitude, rather than power values. Not doing anything to it...\n', cfg.parameter)
        
        if isfield(stimspectra,'powspctrm'), stimspectra = rmfield(stimspectra,'powspctrm'); fprintf('Removing powspctrm field...\n'); end
        if isfield(basespectra,'powspctrm'), basespectra = rmfield(basespectra,'powspctrm'); fprintf('Removing powspctrm field...\n'); end
        
    elseif strcmp(cfg.parameter,'ampspctrm') && isfield(stimspectra,'powspctrm') && isfield(basespectra,'powspctrm')
        fprintf('The field ''%s'' will be calculated by square-rooting the values in ''%s'', which are assumed to be power values.\n', cfg.parameter,'powspctrm')
        for trial = 1:size(stimspectra.powspctrm,1)
            for chan = 1:size(stimspectra.powspctrm,2)
                stimspectra.ampspctrm(trial,chan,:) = sqrt(squeeze(stimspectra.powspctrm(trial,chan,:)));
            end
        end
        for trial = 1:size(basespectra.powspctrm,1)
            for chan = 1:size(basespectra.powspctrm,2)
                basespectra.ampspctrm(trial,chan,:) = sqrt(squeeze(basespectra.powspctrm(trial,chan,:)));
            end
        end
        
        stimspectra = rmfield(stimspectra,'powspctrm'); fprintf('Removing powspctrm field...\n');
        basespectra = rmfield(basespectra,'powspctrm'); fprintf('Removing powspctrm field...\n');
        
    elseif strcmp(cfg.parameter,'ampspctrm')
        error('The cfg option parameter = ''%s'' can only be specified if either ''%s'' or ''%s'' are fields of both input data structures', cfg.parameter,'ampspctrm','powspctrm')
        
    end
    
    varargout{1} = stimspectra;% rmfield(stimspectra,'powspctrm');
    varargout{2} = basespectra;% rmfield(basespectra,'powspctrm');
    varargout{3} = cfg;
    
end

end

