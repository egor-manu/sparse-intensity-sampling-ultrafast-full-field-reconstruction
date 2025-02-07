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

classdef SparseModeDecompositionSimulator < handle
    % SparseModeDecompositionSimulator simulates fiber LP modes, generates random weights,
    % calculates amplitude and intensity patterns, and recovers the weights.

    properties
        modeRelatedParams       % structure of parameters related to the fiber/modes
        dataRelatedParams       % structure of parameters related to the data/image
        LPmodes                 % structure array holding the LP modes
        sensorPositionsMatrix   % matrix that encodes the sensor positions
        Xi                      % matrix containing sparse sampled pairwise mode products
        weightsTrue             % randomly generated mode weights (stored here to plot)
        amplitudePatterns       % field distributions corresponding to weightsTrue
        intensityPatterns       % intensity distributions corresponding to weightsTrue
        weightsRecov            % weights recovered from sparsely taken intensity measurements
    end

    methods
        function obj = SparseModeDecompositionSimulator(modeParams, dataParams)
            obj.modeRelatedParams = modeParams;
            obj.dataRelatedParams = dataParams;

            % find LP modes
            obj.LPmodes = obj.calculateLPmodes(obj.modeRelatedParams, obj.dataRelatedParams);
            assert(~isempty(obj.LPmodes), 'no LP modes found!');
            fprintf('%d LP modes found\n', length(obj.LPmodes));
            
            % Calculate pairwise mode products
            obj.calculatePairwiseModeProducts();

            % Calculate sparse measurement data
            obj.calculateSparseMeasurementsData();
        end

        function plotResults(obj, plotIdx)
            assert(plotIdx>0 && plotIdx<=obj.dataRelatedParams.numberOfSamples, 'sample index must be >0 and < number of samples');
            [~, sensorPositions] = max(obj.sensorPositionsMatrix, [], 2);
            Nx = obj.dataRelatedParams.imageSideLengthPixels;
            Nm = length(obj.LPmodes);
            Cs = zeros(Nx^2, 1);
            Cs(sensorPositions) = 1;
            Cp = obj.reshape2Dto3D(Cs);

            fntSize = 12;

            f = figure('units', 'normalized');
            f.Position = [0.2281 0.3870 0.6490 0.3222];
            tiledlayout('flow','TileSpacing', 'tight', 'Padding', 'none');

            % Plot 1: Noisy intensity and sensor placements
            nexttile;
            pcolor(obj.intensityPatterns(:,:,plotIdx));
            pbaspect([1 1 1]);
            ax1 = gca();
            ax1.XAxis.Visible = 'off';
            ax1.YAxis.Visible = 'off';
            shading flat;
            colormap(ax1, 'bone');
            hold on;
            [rows, cols] = find(Cp == 1);
            plot(cols, rows, 'square', 'LineWidth',3, 'Color',[1,0,1], 'MarkerSize',6);
            title({'Noisy intensity and'; 'sensor placements'}, 'FontSize', fntSize);
            set(gcf, 'Color', [1 1 1]);

            % Plot 2: True intensity (saturation) and phase (color-coded)
            nexttile;
            E = obj.amplitudePatterns(:,:,plotIdx);
            colorData = angle(E);
            alphaData = abs(E).^2/max(abs(E(:)))^2;
            h = surf(zeros(size(E)), colorData, 'EdgeColor', 'none', 'FaceColor', 'interp');
            view(2); pbaspect([1,1,1]); grid off;
            ax2 = gca();
            ax2.XAxis.Visible = 'off';
            ax2.YAxis.Visible = 'off';
            ax2.Color = 0.2*[1 1 1];
            colormap(ax2, 'hsv');
            cl2 = ax2.CLim;
            xlim([1, size(E,1)]);
            ylim([1, size(E,2)]);
            set(h, 'AlphaData', alphaData, 'FaceAlpha', 'interp');
            title({'True intensity (saturation)'; 'and phase (color-coded)'}, 'FontSize', fntSize);

            % Plot 3: Recovered intensity (saturation) and phase (color-coded)
            nexttile;
            amplitudePatternsMeas = obj.calculateAmplitudeForGivenWeights(obj.weightsRecov);
            E = amplitudePatternsMeas(:,:,plotIdx);
            colorData = angle(E);
            alphaData = abs(E).^2/max(abs(E(:)))^2;
            h = surf(zeros(size(E)), colorData, 'EdgeColor', 'none', 'FaceColor', 'interp');
            view(2); pbaspect([1,1,1]); grid off;
            xlim([1, size(E,1)]);
            ylim([1, size(E,2)]);
            clim(cl2);
            ax3 = gca();
            ax3.Color = 0.2*[1 1 1];
            colormap(ax3, 'hsv');
            ax3.XAxis.Visible = 'off';
            ax3.YAxis.Visible = 'off';
            set(h, 'AlphaData', alphaData, 'FaceAlpha', 'interp');
            title({'Retrieved intensity (saturation)'; 'and phase (color-coded)'}, 'FontSize', fntSize);

            % Plot 4: True vs. recovered complex-valued weights
            cWeightsTrue = obj.weightsTrue(plotIdx, 1:Nm) .* exp(1i * [0, obj.weightsTrue(plotIdx, Nm+1:end)]);
            cweightsRecov = obj.weightsRecov(plotIdx, 1:Nm) .* exp(1i * [0, obj.weightsRecov(plotIdx, Nm+1:end)]);
            nexttile;
            plot(cWeightsTrue, 'o', 'LineWidth',2);
            pbaspect([1 1 1]);
            ax = gca();
            ax.YAxis.Direction = 'normal';
            grid on;
            hold on;
            xlim(1.1*[-1, 1]);
            ylim(1.1*[-1, 1]);
            plot(cweightsRecov, 'x', 'LineWidth',2);
            color_o = [0 0.4470 0.7410];
            color_x = [0.8500 0.3250 0.0980];
            ax.FontSize = fntSize;
            ax.YAxis.TickValues = [-1 0 1];
            ax.YAxis.TickLabels = {'-1\it{i}', '0', '1\it{i}'};
            ax.FontWeight = "bold";
            t = title({['\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', color_o) '} True o', ...
                ' \color{black}and \color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', color_x) '} recovered x']; ...
                '\color{black} complex-valued weights'}, ...
                'fontsize', fntSize, 'Interpreter', 'tex');
        end

        function weights = generateRandomWeights(obj)
            Nm = length(obj.LPmodes);
            nSamples = obj.dataRelatedParams.numberOfSamples;
            weights = rand(nSamples, 2*Nm - 1);
            weights = 0.1 + 0.9 * weights;
            weights(:, Nm+1:end) = 2*pi * weights(:, Nm+1:end);
            weights(:, Nm+1) = pi * rand(nSamples, 1); % ensure phi2 is in [0,pi]
            obj.weightsTrue = weights;
        end

        function amplitudePatterns = calculateAmplitudeForGivenWeights(obj, weights)
            Nm = length(obj.LPmodes);
            Nx = obj.dataRelatedParams.imageSideLengthPixels;
            modePatterns = zeros(Nx, Nx, Nm);

            for modeIdx = 1:Nm
                modePatterns(:,:,modeIdx) = obj.LPmodes(modeIdx).Ey;
            end

            weights = reshape(weights, [size(weights, 1), size(weights, 2), 1]); % add dummy dimension (for broadcasting over the 3rd dimension)
            weights = permute(weights, [3, 2, 1]); % swap 1st and 3rd dimensions (for broadcasting over the 3rd dimension)

            amplitudePatterns = modePatterns(:,:,1) .* weights(1,1,:);
            for modeIdx = 2:Nm
                cWeights = weights(1, modeIdx, :).*exp(1i*weights(1, Nm+modeIdx-1, :));
                amplitudePatterns = amplitudePatterns + modePatterns(:,:,modeIdx) .* cWeights;
            end
        end

        function intensityPatterns = calculateIntensityForGivenWeights(obj, weights)
            obj.amplitudePatterns = calculateAmplitudeForGivenWeights(obj, weights);
            intensityPatterns = abs(obj.amplitudePatterns).^2;

            % Add noise to the intensity patterns
            maxIntensity = max(intensityPatterns(:));
            intensityPatterns = intensityPatterns + maxIntensity * obj.dataRelatedParams.relativeNoise * randn(size(intensityPatterns));
            obj.intensityPatterns = intensityPatterns;
        end

        function array2D = reshape3Dto2D(~, array3D)
            [dim1, ~, dim3] = size(array3D);
            array2D = reshape(array3D, [dim1 * dim1, dim3]);
        end

        function weightsRecov = recoverWeights(obj, weightScalarProducts)
            Nm = length(obj.LPmodes);
            nSamples = size(weightScalarProducts, 2);
            weightsRecov = zeros(nSamples, 2*Nm - 1);
            weightsRecov(:, 1:Nm) = sqrt(abs(weightScalarProducts(1:Nm, :))).';
            for idx = 1:Nm-1
                cosVals = weightScalarProducts(Nm+idx, :)' ./ weightsRecov(:,1) ./ weightsRecov(:, idx+1) / 2;
                cosVals(cosVals > 1) = 1;
                cosVals(cosVals < -1) = -1;
                weightsRecov(:, Nm+idx) = acos(cosVals);
            end
            for idx = 3:Nm
                cosOther = cos(weightsRecov(:, Nm+1) + weightsRecov(:, Nm-1+idx));
                cosCurrent = cos(weightsRecov(:, Nm+1) - weightsRecov(:, Nm-1+idx));
                trueVals = weightScalarProducts(2*Nm-3+idx, :)' ./ weightsRecov(:,2) ./ weightsRecov(:, idx) / 2;
                signRemains = -2*(abs(cosOther - trueVals) < abs(cosCurrent - trueVals)) + 1;
                weightsRecov(:, Nm-1+idx) = mod(weightsRecov(:, Nm-1+idx) .* signRemains, 2*pi);
            end
            obj.weightsRecov = weightsRecov;
        end

    end

    methods(Access = private)
        function LPmodes = calculateLPmodes(~, modeRelatedParams, dataRelatedParams)
            Nx = dataRelatedParams.imageSideLengthPixels;
            fiberDtoImgWidth = dataRelatedParams.fiberDtoImgWidthRatio;
            fiberCoreRadius = 1e-6 * modeRelatedParams.fiberCoreDiameter_um / 2;
            wavelength = 1e-9 * modeRelatedParams.wavelength_nm;
            NA = modeRelatedParams.numericalAperture;
            n = modeRelatedParams.claddingRefractiveIndex;
            R = fiberCoreRadius / fiberDtoImgWidth;
            dn = NA^2 / (2*n);
            nCore = n + dn;
            nClad = n;
            k0 = 2*pi/wavelength;
            v = (k0 * fiberCoreRadius) * sqrt(nCore^2 - nClad^2);
            us = linspace(0, 0.999*v, 1e4);
            L = 0;
            x = linspace(-R, R, Nx);
            y = linspace(-R, R, Nx);
            [X, Y] = ndgrid(x, y);
            [PHI, RO] = cart2pol(X, Y);

            % Preallocate structure array for LP modes
            LPmodes = struct('name', {}, 'name_plot', {}, 'L', {}, 'M', {}, ...
                'betar', {}, 'u', {}, 'w', {}, 'Ey', {}, 'x', {}, 'y', {});

            modeNumber = 1;
            while true
                dispRelFunc = @(x) x/nCore .* besselj(L-1, x) .* besselk(L, sqrt(v^2 - x.^2)) + ...
                    sqrt(v^2 - x.^2)/nClad .* besselk(L-1, sqrt(v^2 - x.^2)) .* besselj(L, x);
                f = dispRelFunc(us);
                z = f .* circshift(f, [0, -1]);
                z(end) = z(end-1);
                uModes = find(z < 0);
                if ~nnz(z < 0)
                    break
                end
                for M = 1:length(uModes)
                    try
                        u = fzero(dispRelFunc, [us(uModes(M)-1), us(uModes(M)+1)]);
                    catch ex
                        disp(ex.message)
                        continue
                    end
                    w = sqrt(v^2 - u^2);
                    betar = sqrt(k0^2 * fiberCoreRadius^2 * nCore^2 - u^2);
                    currentMode.L = L;
                    currentMode.M = M;
                    currentMode.u = u;
                    currentMode.w = w;
                    currentMode.betar = betar;
                    currentMode.x = x;
                    currentMode.y = y;

                    Ey = zeros(size(X));
                    in_core = RO <= fiberCoreRadius;
                    in_clad = RO > fiberCoreRadius;
                    Ey_bessel_in = besselj(L, u * RO/fiberCoreRadius) / besselj(L, u);
                    Ey_bessel_out = besselk(L, w * RO/fiberCoreRadius) / besselk(L, w);

                    cosPart = cos(L * PHI);
                    Ey(in_core) = Ey_bessel_in(in_core) .* cosPart(in_core);
                    Ey(in_clad) = Ey_bessel_out(in_clad) .* cosPart(in_clad);
                    currentMode.Ey = Ey / max(abs(Ey(:)));
                    currentMode.name = ['LP', num2str(L), num2str(M)];
                    currentMode.name_plot = ['LP_{', num2str(L), num2str(M), '}'];

                    if L == 0
                        LPmodes(modeNumber) = currentMode;
                        modeNumber = modeNumber + 1;
                    else
                        % For nonzero L, generate cosine and sine versions
                        currentMode.name = [currentMode.name, 'c'];
                        currentMode.name_plot = [currentMode.name_plot, '^c'];
                        LPmodes(modeNumber) = currentMode;
                        modeNumber = modeNumber + 1;

                        sinPart = sin(L * PHI);
                        Ey(in_core) = Ey_bessel_in(in_core) .* sinPart(in_core);
                        Ey(in_clad) = Ey_bessel_out(in_clad) .* sinPart(in_clad);

                        currentMode.Ey = Ey / max(abs(Ey(:)));
                        currentMode.name = ['LP', num2str(L), num2str(M), 's'];
                        currentMode.name_plot = ['LP_{', num2str(L), num2str(M), '}^s'];
                        LPmodes(modeNumber) = currentMode;
                        modeNumber = modeNumber + 1;
                    end
                end
                L = L + 1;
            end
        end

        function array3D = reshape2Dto3D(~, array2D)
            [numPix, dim3] = size(array2D, 1);
            dim1 = sqrt(numPix);
            array3D = reshape(array2D, [dim1, dim1, dim3]);
        end

        function pairwiseModeProducts = calculatePairwiseModeProducts(obj)
            Nx = obj.dataRelatedParams.imageSideLengthPixels;
            Nm = length(obj.LPmodes);
            vectorizedModes = zeros(Nx*Nx, Nm);
            pairwiseModeProducts = zeros(Nx*Nx, Nm*(Nm+1)/2);
            for modeIdx = 1:Nm
                vectorizedModes(:, modeIdx) = reshape(obj.LPmodes(modeIdx).Ey, [Nx*Nx, 1]);
                pairwiseModeProducts(:, modeIdx) = vectorizedModes(:, modeIdx).^2;
            end
            idxk = Nm;
            for modeIdx1 = 1:Nm
                for modeIdx2 = modeIdx1+1:Nm
                    idxk = idxk + 1;
                    pairwiseModeProducts(:, idxk) = vectorizedModes(:, modeIdx1) .* vectorizedModes(:, modeIdx2);
                end
            end
        end

        function calculateSparseMeasurementsData(obj)
            pairwiseModeProducts = obj.calculatePairwiseModeProducts();
            [pwpOrth, pwpOrthTRransform] = qr(pairwiseModeProducts, 'econ');
            nSensors = size(pairwiseModeProducts, 2);
            [~, ~, pivots] = qr(pwpOrth', 'vector', 'econ');
            obj.sensorPositionsMatrix = zeros(nSensors, size(pairwiseModeProducts, 1));
            for sensorIdx = 1:nSensors
                obj.sensorPositionsMatrix(sensorIdx, pivots(sensorIdx)) = 1;
            end
            Theta = obj.sensorPositionsMatrix * pwpOrth;
            obj.Xi = pinv(pwpOrthTRransform) * pinv(Theta);
        end
    end
end
