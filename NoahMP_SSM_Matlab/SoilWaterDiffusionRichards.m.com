function [MatRight, MatLeft1, MatLeft2, MatLeft3, DrainSoilBot, SoilWatConductivity, SoilWatDiffusivity] = SoilWaterDiffusionRichards(noahmp)
% SoilWaterDiffusionRichards: Solve Richards equation for soil water movement/diffusion
% Compute the right hand side of the time tendency term of the soil
% water diffusion equation. Also prepares the matrix
% coefficients for the tri-diagonal matrix of the implicit time scheme.
% This is simplified for Runoff Option 3 and OptPermeability1
% Initialize output variables

NumSoilLayer = noahmp.config.domain.NumSoilLayer;
MatRight = zeros(NumSoilLayer, 1);
MatLeft1 = zeros(NumSoilLayer, 1);
MatLeft2 = zeros(NumSoilLayer, 1);
MatLeft3 = zeros(NumSoilLayer, 1);
SoilWatConductivity = zeros(NumSoilLayer, 1);
SoilWatDiffusivity = zeros(NumSoilLayer, 1);

% Initialize local variables
DepthSnowSoilInv = zeros(NumSoilLayer, 1);
SoilThickTmp = zeros(NumSoilLayer, 1);
SoilWaterGrad = zeros(NumSoilLayer, 1);
WaterExcess = zeros(NumSoilLayer, 1);
SoilMoistureTmp = zeros(NumSoilLayer, 1);

% Compute soil hydraulic conductivity and diffusivity
for LoopInd = 1:NumSoilLayer
    [SoilWatDiffusivity(LoopInd), SoilWatConductivity(LoopInd)] = SoilDiffusivityConductivity(noahmp, ...
        noahmp.water.state.SoilMoisture(LoopInd), noahmp.water.state.SoilImpervFrac(LoopInd), LoopInd,1);
    SoilMoistureTmp(LoopInd) = noahmp.water.state.SoilMoisture(LoopInd);
end

% Compute gradient and flux of soil water diffusion terms
for LoopInd = 1:NumSoilLayer
    if LoopInd == 1
        SoilThickTmp(LoopInd) = -noahmp.config.domain.DepthSoilLayer(LoopInd);
        DepthSnowSoilTmp = -noahmp.config.domain.DepthSoilLayer(LoopInd+1);
        DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp;
        SoilWaterGrad(LoopInd) = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp;
        WaterExcess(LoopInd) = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - ...
            noahmp.water.flux.InfilRateSfc + noahmp.water.flux.TranspWatLossSoilMean(LoopInd) + noahmp.water.flux.EvapSoilSfcLiqMean;
    elseif LoopInd < NumSoilLayer
        SoilThickTmp(LoopInd) = noahmp.config.domain.DepthSoilLayer(LoopInd-1) - noahmp.config.domain.DepthSoilLayer(LoopInd);
        DepthSnowSoilTmp = noahmp.config.domain.DepthSoilLayer(LoopInd-1) - noahmp.config.domain.DepthSoilLayer(LoopInd+1);
        DepthSnowSoilInv(LoopInd) = 2.0 / DepthSnowSoilTmp;
        SoilWaterGrad(LoopInd) = 2.0 * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd+1)) / DepthSnowSoilTmp;
        WaterExcess(LoopInd) = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) - ...
            SoilWatDiffusivity(LoopInd-1) * SoilWaterGrad(LoopInd-1) - SoilWatConductivity(LoopInd-1) + ...
            noahmp.water.flux.TranspWatLossSoilMean(LoopInd);
    else
        SoilThickTmp(LoopInd) = noahmp.config.domain.DepthSoilLayer(LoopInd-1) - noahmp.config.domain.DepthSoilLayer(LoopInd);
        DrainSoilBot = noahmp.water.param.SoilDrainSlope * SoilWatConductivity(LoopInd);
        WaterExcess(LoopInd) = -(SoilWatDiffusivity(LoopInd-1) * SoilWaterGrad(LoopInd-1)) - SoilWatConductivity(LoopInd-1) + ...
            noahmp.water.flux.TranspWatLossSoilMean(LoopInd) + DrainSoilBot;
    end
end

% Prepare the matrix coefficients for the tri-diagonal matrix
for LoopInd = 1:NumSoilLayer
    if LoopInd == 1
        MatLeft1(LoopInd) = 0.0; % A In thr Matrix Solver
        MatLeft2(LoopInd) = SoilWatDiffusivity(LoopInd) * DepthSnowSoilInv(LoopInd) / SoilThickTmp(LoopInd); % B in the Matrix Solver
        MatLeft3(LoopInd) = -MatLeft2(LoopInd); % C in the Matrix Solver
    elseif LoopInd < NumSoilLayer
        MatLeft1(LoopInd) = -SoilWatDiffusivity(LoopInd-1) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd);
        MatLeft3(LoopInd) = -SoilWatDiffusivity(LoopInd) * DepthSnowSoilInv(LoopInd) / SoilThickTmp(LoopInd);
        MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    else
        MatLeft1(LoopInd) = -SoilWatDiffusivity(LoopInd-1) * DepthSnowSoilInv(LoopInd-1) / SoilThickTmp(LoopInd);
        MatLeft3(LoopInd) = 0.0;
        MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    end
    MatRight(LoopInd) = WaterExcess(LoopInd) / (-SoilThickTmp(LoopInd)); % D in the Matrix Solver
end

end
