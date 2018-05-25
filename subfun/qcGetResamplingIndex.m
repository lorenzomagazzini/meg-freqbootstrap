function [ bootindex, bootindexstim, bootindexbase ] = qcGetResamplingIndex( numboot, numtrial, numtrialstim, numtrialbase )
%[ bootindex, bootindexstim, bootindexbase ] = qcGetResamplingIndex( numboot, numtrial, numtrialstim, numtrialbase )
%   Detailed explanation goes here

%one condition
if numtrial>0 && ~istrue(any([numtrialstim, numtrialbase]))
    
    %prepare the array to store the resampling indices
    bootindex = nan(numboot,numtrial);
    bootindexstim = [];
    bootindexbase = [];
    
    %get the resampling indices
    for boot = 1:numboot
        bootindex(boot,:) = randsample(numtrial, numtrial, 'true');
    end
    
%two conditions
elseif ~istrue(numtrial) && numtrialstim>0 && numtrialbase>0
    
    %prepare the array to store the resampling index
    bootindex = [];
    bootindexstim = nan(numboot,numtrialstim);
    bootindexbase = nan(numboot,numtrialbase);
    
    %same number of trials for stimulus and baseline
    if numtrialstim == numtrialbase
        for boot = 1:numboot
            %get the same resampling indices
            bootindexstim(boot,:) = randsample(numtrialstim, numtrialstim, 'true');
            bootindexbase = bootindexstim;
        end
        %different number of trials for stimulus and baseline
    elseif numtrialstim ~= numtrialbase
        for boot = 1:numboot
            %get different resampling indices
            bootindexstim(boot,:) = randsample(numtrialstim, numtrialstim, 'true');
            bootindexbase(boot,:) = randsample(numtrialbase, numtrialbase, 'true');
        end
    end
    
end

end

