function CtrCommand = StateControl(x_Nom, x_Ref, u_Ref, ...
                        x_Ctrl, iTime, CtrlStruct, IdStruct, SimParam)
    dt    = 1/SimParam.SamplingFreq;
    IdXRW = IdStruct.RealMdl.IdXRW;
    
    CtrCommand = u_Ref(:,iTime)-CtrlStruct.Kd*(x_Ctrl-(x_Nom+x_Ref(:,iTime)));
    % max torque 
    CtrCommand = max(min(CtrCommand,...,inf),-inf);...
                        SimParam.MaxTorqueRW),...
                        -SimParam.MaxTorqueRW);
    % clamp RW speeds CtrlStruct.InertiaRW 
    CtrCommand = max(min(CtrCommand,...,inf),-inf);...
             CtrlStruct.InertiaRW.*(CtrlStruct.maxrwspeed - x_Ctrl(IdXRW))/dt),...
             CtrlStruct.InertiaRW.*(CtrlStruct.minrwspeed - x_Ctrl(IdXRW))/dt);
end

