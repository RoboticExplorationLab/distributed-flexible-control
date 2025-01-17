function [OBSWStruct, ModelStruct, IdStruct] = LQRDesigner(SysIddModel, ...
                                OBSWStruct, ModelStruct, IdStruct, SimParam)
    % Infinite Horizon LQR
    %% Load Vars
    Nmodr        = IdStruct.RealMdl.Nmod;
    RealID       = IdStruct.RealMdl;
    IdKD         = IdStruct.IdKD;
    Ad = SysIddModel.A;
    Bd = SysIddModel.B(:,IdStruct.IdUWU);
    CtrlStruct              = struct;
    NX               = RealID.NFilt + RealID.Nrig + 2*Nmodr ...
                                                            + RealID.NXRW;
    CtrlStruct.IdXRW = RealID.NFilt + RealID.Nrig ... 
                                    + 2*Nmodr + (1:RealID.NXRW);
    NU                 = IdStruct.NU;  
    CtrlStruct.NU      = NU;
    CtrlStruct.NX      = NX;
  
    Q  = SimParam.QCtrl;
    R  = SimParam.RCtrl;

    %% Infinite LQR
    nxctrl = IdKD.NCOaOB;
    T2Ctrl = OBSWStruct.TKD(:,IdKD.IdCOaOB)';
    Tact = OBSWStruct.TKD(IdStruct.RealMdl.IdXRW,[IdKD.IdXMomRW, IdKD.IdXNullRWd])';
    Ad = T2Ctrl * Ad * T2Ctrl';
    Bd = T2Ctrl * Bd * Tact';
    [Kd1,P,CtrlPoles] = dlqr(Ad, Bd, Q , R, zeros(nxctrl, NU));
    Kd    = Tact' * Kd1 * T2Ctrl;
    
    %% Save Vars
    OBSWStruct.Ctrl           = CtrlStruct;
    OBSWStruct.Ctrl.Kd        = Kd;
    OBSWStruct.Ctrl.P         = P;
    OBSWStruct.Ctrl.T2Ctrl    = T2Ctrl;
    OBSWStruct.Ctrl.Tact      = Tact; 
    OBSWStruct.Ctrl.CtrlPoles = CtrlPoles;
    OBSWStruct.Ctrl.Ad        = Ad;
    OBSWStruct.Ctrl.Bd        = Bd;
    OBSWStruct.Ctrl.InertiaRW = ModelStruct.Act.InertiaRW;
    OBSWStruct.Ctrl.minrwspeed = ModelStruct.Act.minrwspeed;
    OBSWStruct.Ctrl.nomrwspeed = ModelStruct.Act.nomrwspeed;
    OBSWStruct.Ctrl.maxrwspeed = ModelStruct.Act.maxrwspeed;
end
