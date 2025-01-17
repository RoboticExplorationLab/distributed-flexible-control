function [x_ref, u_ref] = StateGuidance(tVec, RealSysCont,...
                                        SimParam, OBSWStruct, IdStruct)

    if strcmp(SimParam.GuidMethod,'ZVTOrig')
        [x_ref, u_ref] = ZVTORigSlewGuidance( tVec, RealSysCont, SimParam, ...
                                                        OBSWStruct, IdStruct);
    else
        x_ref = zeros(IdStruct.RealMdl.NX, length(tVec));
        u_ref = zeros(IdStruct.NU, length(tVec));
    end
end
