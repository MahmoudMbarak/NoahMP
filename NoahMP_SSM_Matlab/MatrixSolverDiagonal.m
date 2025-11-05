function [P, Delta] = MatrixSolverTriDiagonal(A, B, C, D, IndTopLayer, NumSoilLayer, NumSnowLayerMax)
% MatrixSolverTriDiagonal: Solve tri-diagonal matrix problem
% ------------------------ Code history --------------------------------------------------
% Original Noah-MP subroutine: ROSR12
% Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
% Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
% Translated to MATLAB: [Mahmoud Mbarak]
% ----------------------------------------------------------------------------------------

% Initialize P and Delta to match the size of input arrays
P = zeros(size(A));
Delta = zeros(size(A));

% Initialize EQN COEF C for the lowest soil layer
C(NumSoilLayer) = 0.0;
P(IndTopLayer) = -C(IndTopLayer) / B(IndTopLayer);

% Solve the coefs for the 1st soil layer
Delta(IndTopLayer) = D(IndTopLayer) / B(IndTopLayer);

% Solve the coefs for soil layers 2 thru NumSoilLayer
for K = IndTopLayer+1:NumSoilLayer
    P(K) = -C(K) * (1.0 / (B(K) + A(K) * P(K-1)));
    Delta(K) = (D(K) - A(K) * Delta(K-1)) * (1.0 / (B(K) + A(K) * P(K-1)));
end

% Set P to Delta for lowest soil layer
P(NumSoilLayer) = Delta(NumSoilLayer);

% Adjust P for soil layers 2 thru NumSoilLayer
for K = IndTopLayer+1:NumSoilLayer
    KK = NumSoilLayer - K + (IndTopLayer-1) + 1;
    P(KK) = P(KK) * P(KK+1) + Delta(KK);
end

end
