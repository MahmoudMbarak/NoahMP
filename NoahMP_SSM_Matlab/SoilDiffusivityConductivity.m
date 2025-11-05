function [SoilWatDiffusivity, SoilWatConductivity] = SoilDiffusivityConductivity(noahmp, SoilMoisture, SoilImpervFrac, IndLayer, option)
% SoilDiffusivityConductivityOpt1: Calculate soil water diffusivity and hydraulic conductivity based on Clapp and Hornberger
% This function implements Option 1: linear effects (more permeable, Niu and Yang, 2006)

if option == 1
    % Extract necessary parameters from noahmp structure
    SoilMoistureSat = noahmp.water.param.SoilMoistureSat(IndLayer);
    SoilExpCoeffB = noahmp.water.param.SoilExpCoeffB(IndLayer);
    SoilWatDiffusivitySat = noahmp.water.param.SoilWatDiffusivitySat(IndLayer);
    SoilWatConductivitySat = noahmp.water.param.SoilWatConductivitySat(IndLayer);
    % Calculate soil pre-factor
    SoilPreFac = max(0.01, SoilMoisture / SoilMoistureSat);

    % Calculate soil water diffusivity
    SoilExpTmp = SoilExpCoeffB + 2.0;
    SoilWatDiffusivity = SoilWatDiffusivitySat * SoilPreFac ^ SoilExpTmp;
    SoilWatDiffusivity = SoilWatDiffusivity * (1.0 - SoilImpervFrac);

    % Calculate soil hydraulic conductivity
    SoilExpTmp = 2.0 * SoilExpCoeffB + 3.0;
    SoilWatConductivity = SoilWatConductivitySat * SoilPreFac ^ SoilExpTmp;
    SoilWatConductivity = SoilWatConductivity * (1.0 - SoilImpervFrac);

else
    %% Replicating the Implicit scheme:
    SoilConductivityFactA = 1.175 * 10^ 6;
    SoilConductivityFactAlpha = 1.611 * 10^6;
    SoilConductivityFactBeta = 3.96;
    
    SoilMoistureSat = noahmp.water.param.SoilMoistureSat(IndLayer);
    SoilMoistureResidual = 0.075;
    SoilWatConductivitySat = noahmp.water.param.SoilWatConductivitySat(IndLayer);


    SoilConductivityFactGamma = 4.74; 
    SoilPsi = (SoilConductivityFactAlpha .* abs(SoilMoistureSat - SoilMoisture) ./ abs(SoilMoisture - SoilMoistureResidual)).^(1/SoilConductivityFactBeta);
    SoilWatConductivity = (SoilWatConductivitySat*SoilConductivityFactA)./(SoilConductivityFactA+(abs(SoilPsi).^SoilConductivityFactGamma));

    SoildPsidSoilMoisture = (1./SoilConductivityFactBeta) .* SoilConductivityFactAlpha.*(SoilMoistureSat - SoilMoistureResidual) ...
        .*(SoilConductivityFactAlpha.*abs(SoilMoistureSat - SoilMoisture) ...
        ./abs(SoilMoisture - SoilMoistureResidual)) ...
        .^(1/SoilConductivityFactBeta - 1) ./ ...
        (SoilMoisture - SoilMoistureResidual).^2;

    SoilWatDiffusivity = SoildPsidSoilMoisture * SoilWatConductivity;

end
