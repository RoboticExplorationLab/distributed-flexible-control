function [wKnown,wUnknown] = RWPertFunc(thetarw,wrw,ModelStruct,IdStruct,SimParam)
        NWRWF             = IdStruct.NWRWM;
        NWRWM             = IdStruct.NWRWF;
        IdUWWRWF          = IdStruct.IdUWWRWF;
        IdUWWRWM          = IdStruct.IdUWWRWM;
        IdUWWRW           = IdStruct.IdUWWRW;
        NWPerNU           = 2*SimParam.RadFPertRW+2*SimParam.RadMPertRW;
        NU                = IdStruct.NU;
        PertMag           = ones(NWRWM+NWRWF,1);
        PertMag(IdUWWRWF) = ModelStruct.Pert.FRWRadFHCoeff;
        PertMag(IdUWWRWM) = ModelStruct.Pert.MRWRadFHCoeff;
        PertMag      = PertMag .* repelem(wrw.^ 2,NWPerNU,1);
        PertPhas     = repmat(deg2rad([0,90,-90,0]),1,NU)';
        thetarwpert  = repelem(thetarw,NWPerNU,1);

        wKnown(IdUWWRW) = PertMag .* sin(thetarwpert + PertPhas);
        wUnknown(IdUWWRWF) = SimParam.FRWRadBBVar.^0.5*randn(NWRWF,1);
        wUnknown(IdUWWRWM) = SimParam.MRWRadBBVar.^0.5*randn(NWRWM,1);    
end

