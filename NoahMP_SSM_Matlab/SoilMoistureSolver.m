function noahmp = SoilMoistureSolver(noahmp, TimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)
% SoilMoistureSolver: Compute soil moisture content based on Richards diffusion & tri-diagonal matrix

% Extract necessary variables from noahmp structure
NumSoilLayer = noahmp.config.domain.NumSoilLayer;
DepthSoilLayer = noahmp.config.domain.DepthSoilLayer;
ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
SoilMoistureSat = noahmp.water.param.SoilMoistureSat;
SoilIce = noahmp.water.state.SoilIce;
SoilLiqWater = noahmp.water.state.SoilLiqWater;
SoilMoisture = noahmp.water.state.SoilMoisture;

% Initialize local variables
SoilSaturationExcess = 0;
SoilEffPorosity = zeros(NumSoilLayer, 1);

% Update tri-diagonal matrix elements
MatRight = MatRight * TimeStep;
MatLeft1 = MatLeft1 * TimeStep;
MatLeft2 = 1.0 + MatLeft2 * TimeStep;
MatLeft3 = MatLeft3 * TimeStep;

% Call MatrixSolverTriDiagonal to solve the tri-diagonal matrix
[DeltaTheta, ~] = MatrixSolverTriDiagonal(MatLeft1, MatLeft2, MatLeft3, MatRight, 1, NumSoilLayer, 0);

% Update SoilLiqWater
SoilLiqWater = SoilLiqWater + DeltaTheta;

% Excessive water above saturation in a layer is moved to its unsaturated layer like in a bucket
for LoopInd = NumSoilLayer:-1:2
    SoilEffPorosity(LoopInd) = max(1.0e-4, (SoilMoistureSat(LoopInd) - SoilIce(LoopInd)));
    SoilSaturationExcess = max((SoilLiqWater(LoopInd) - SoilEffPorosity(LoopInd)), 0) * ThicknessSnowSoilLayer(LoopInd);
    SoilLiqWater(LoopInd) = min(SoilEffPorosity(LoopInd), SoilLiqWater(LoopInd));
    SoilLiqWater(LoopInd-1) = SoilLiqWater(LoopInd-1) + SoilSaturationExcess / ThicknessSnowSoilLayer(LoopInd-1);
end

SoilEffPorosity(1) = max(1.0e-4, (SoilMoistureSat(1) - SoilIce(1)));
SoilSaturationExcess = max((SoilLiqWater(1) - SoilEffPorosity(1)), 0) * ThicknessSnowSoilLayer(1);
SoilLiqWater(1) = min(SoilEffPorosity(1), SoilLiqWater(1));

if SoilSaturationExcess > 0
    SoilLiqWater(2) = SoilLiqWater(2) + SoilSaturationExcess / ThicknessSnowSoilLayer(2);
    for LoopInd = 2:NumSoilLayer-1
        SoilEffPorosity(LoopInd) = max(1.0e-4, (SoilMoistureSat(LoopInd) - SoilIce(LoopInd)));
        SoilSaturationExcess = max((SoilLiqWater(LoopInd) - SoilEffPorosity(LoopInd)), 0) * ThicknessSnowSoilLayer(LoopInd);
        SoilLiqWater(LoopInd) = min(SoilEffPorosity(LoopInd), SoilLiqWater(LoopInd));
        SoilLiqWater(LoopInd+1) = SoilLiqWater(LoopInd+1) + SoilSaturationExcess / ThicknessSnowSoilLayer(LoopInd+1);
    end
    SoilEffPorosity(NumSoilLayer) = max(1.0e-4, (SoilMoistureSat(NumSoilLayer) - SoilIce(NumSoilLayer)));
    SoilSaturationExcess = max((SoilLiqWater(NumSoilLayer) - SoilEffPorosity(NumSoilLayer)), 0) * ThicknessSnowSoilLayer(NumSoilLayer);
    SoilLiqWater(NumSoilLayer) = min(SoilEffPorosity(NumSoilLayer), SoilLiqWater(NumSoilLayer));
end

% Update SoilMoisture
SoilMoisture = SoilLiqWater + SoilIce;

% Update noahmp structure with new values
noahmp.water.state.SoilLiqWater = SoilLiqWater;
noahmp.water.state.SoilMoisture = SoilMoisture;
noahmp.water.state.SoilEffPorosity = SoilEffPorosity;
noahmp.water.state.SoilSaturationExcess = SoilSaturationExcess;
end
