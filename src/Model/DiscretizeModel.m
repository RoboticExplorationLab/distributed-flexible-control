function [SysMeasRedDisc, IdStruct, ModelStruct] = ...
              DiscretizeModel(SysMeasRed, ModelStruct, SimParam,IdStruct)
    %% Anti-aliasing filter on analog measurement signal
    % sampling frequency
    SamplingFreq = SimParam.SamplingFreq; % [Hz]
    % sampling time
    SamplingTime = 1/SamplingFreq; % [s]
    RealId       = IdStruct.RealMdl;

    %% Discretization of Reduced Order System
    BigA = [SysMeasRed.E\SysMeasRed.A, SysMeasRed.E\SysMeasRed.B;...
            zeros(IdStruct.NUW,RealId.NX+IdStruct.NUW)];
    AA = expm(BigA*SamplingTime);
    A = AA(1:RealId.NX,1:RealId.NX);
    B = AA(1:RealId.NX,RealId.NX+1:end);
    SysMeasRedDisc = ss(A, B, SysMeasRed.C, SysMeasRed.D, SamplingTime);

    ModelStruct.Pert.GDistJacDisc = SysMeasRedDisc.B(:,IdStruct.IdUWW);
end

