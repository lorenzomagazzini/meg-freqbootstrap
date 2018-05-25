function [ bootspectra ] = qcGetResampledSingleTrialSpectra( origspectra, bootindex, parameter )
%[ bootspectra ] = qcGetResampledSingleTrialSpectra( origspectra, bootindex, parameter )
%   See inside function

%simulate fieldtrip structure
bootspectra = origspectra;

%resample single-trial spectra
bootpowspctrm = origspectra.(parameter)(bootindex,:,:);

%replace powspctrm field in ft structure
bootspectra.(parameter) = bootpowspctrm;

end

