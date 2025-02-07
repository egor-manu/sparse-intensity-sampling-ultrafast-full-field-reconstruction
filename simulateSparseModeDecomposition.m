%% ==================== Code Attribution and Citation =====================
%
% This code was developed by Egor Manuylovich and is associated with the research
% described in the following publication:
%
%    Sparse intensity sampling for ultrafast full-field reconstruction in low-dimensional photonic systems
%    DOI: doi.org/10.21203/rs.3.rs-3921498/v1
%
% If you use or adapt any part of this code in your work, please acknowledge its
% source by citing the above publication.
%
% Example citation:
%    E Manuylovich, "Sparse intensity sampling for ultrafast full-field reconstruction in low-dimensional photonic systems". DOI: doi.org/10.21203/rs.3.rs-3921498/v1
%
% Thank you for your interest in and use of this code!
% =========================================================================

clc
clear
close all
rng(41)

% Define mode-related parameters
modeRelatedParams.fiberCoreDiameter_um = 5;         % [um] fiber core diameter
modeRelatedParams.wavelength_nm = 650;              % [nm] wavelength
modeRelatedParams.numericalAperture = 0.14;         % fiber NA
modeRelatedParams.claddingRefractiveIndex = 1.46;   % fiber clad refractive index

% Define data-related parameters
dataRelatedParams.imageSideLengthPixels = 64;       % [pix] size of the simulated intensity patterns
dataRelatedParams.fiberDtoImgWidthRatio = 0.75;     % relative size of fiber core with respect to intensity pattern
dataRelatedParams.numberOfSamples = 10;             % number of simulated intensity patterns
dataRelatedParams.relativeNoise = 1e-2;             % relative intensity noise

% create an instance of the simulation class
s = SparseModeDecompositionSimulator(modeRelatedParams, dataRelatedParams);

% matrix containing sparse sampled pairwise mode products (Xi matrix in the paper)
Xi = s.Xi;

% matrix that encodes the sensor positions (P matrix in the paper)
sensorPositionsMatrix = s.sensorPositionsMatrix;

% generate random mode weights (C in the paper)
weightsTrue = s.generateRandomWeights();

% calculate (noisy) intensity patterns for given weights (I in the paper)
intensityPatterns = s.calculateIntensityForGivenWeights(weightsTrue);   

% pick intensity measurements only from the sensor positions (Ip in the paper)
sparselySampledIntensity = sensorPositionsMatrix * s.reshape3Dto2D(intensityPatterns);

% calculate "intensity weights" (S weights in the paper)
weightScalarProducts = s.Xi * sparselySampledIntensity;

% recover mode weights from "intensity weights" (C in the paper)
weightsRecov = s.recoverWeights(weightScalarProducts);

% calculate MAE
mae = mean(abs(weightsTrue - weightsRecov), 'all');
fprintf('MAE = %1.2e, averaged over %d samples and %d mode weights.\n', mae, size(weightsTrue, 1), size(weightsTrue, 2));

% plot results
sampleIndex2plot = 1;
s.plotResults(sampleIndex2plot)

