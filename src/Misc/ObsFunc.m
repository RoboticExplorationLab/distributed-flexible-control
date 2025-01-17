function [meas,C,D,meas_noerr] = ObsFunc(RealStruct, State, Input, wPert, ...
                                   thetarw, IdStruct, ModelStruct, SimParam)
    % Observation function: receives the model and states, outputs the
    % measurement
    % Taking matrices from ss is slow for some reason
    C          = RealStruct.RealSysDiscC;
    D          = RealStruct.RealSysDiscD;
    meas_noerr = C * State + D * [Input; wPert];
    
    if ~strcmp(SimParam.obsmet,'linear') && IdStruct.NYAMU && ~IdStruct.RealMdl.NFilt
        % nonlinear acceleration terms
        omega  = State(IdStruct.RealMdl.IdXrigrd);
        eta    = State(IdStruct.RealMdl.IdXmodq);
        etadot = State(IdStruct.RealMdl.IdXmodd);
        AMUDir = ModelStruct.Sens.AMUDir';
        ramu0  = ModelStruct.Sens.AMUPos';
        [~,pos]=intersect(ModelStruct.IdInner, ModelStruct.ActNodesID);
        NPhi = zeros(3,IdStruct.NYAMU,IdStruct.RealMdl.Nmod);
        IdUFlex = SimParam.RWperpoint(1)+1;
        if IdStruct.Nctrlpoints > 1
            NPhi(3,IdUFlex:end,:) = ModelStruct.Node2Mode(pos*3-2,:);
        end
        Nq    = tensorprod(NPhi,eta,3,1);
        Nqdot = tensorprod(NPhi,etadot,3,1);
        ramu  = ramu0 + Nq;
        
        RealStruct.FiltFlag = false;
        DynFunc = @(x,u) rigflexcontDisc(x,u, RealStruct, IdStruct,...
                                                    ModelStruct, SimParam);
        NLXA = DynFunc([State;thetarw],Input); % nonlinear state acceleration
        % trans. from state accel. to AMU accel.
        idxamu = [IdStruct.RealMdl.IdXrigd, IdStruct.RealMdl.IdXmodd];
        TXAMU = zeros(IdStruct.NYAMU, length(idxamu));
        for i = 1:IdStruct.NYAMU
            TXAMU(i,:) = AMUDir(:,i)'*[eye(3), -skew(ramu(:,i)), squeeze(NPhi(:,i,:))];
        end
        
        meas_noerr(IdStruct.IdYAMU) = TXAMU * NLXA(idxamu);
        % cent. and coriolis accel on every accelerometer position
        anonlin = skew(omega) * Nqdot + skew(omega) * (skew(omega) * ramu);
        meas_noerr(IdStruct.IdYAMU) = meas_noerr(IdStruct.IdYAMU) ...
                                                   + dot(AMUDir,anonlin,1)';
    end
    
    % % Error model
    meas       = meas_noerr;
    if SimParam.measerr
        if IdStruct.NYAMU
            AMUWNVar    = ModelStruct.Sens.GyroWNVar;
            AMUBias     = ModelStruct.Sens.AMUBias;
            AMUWN       = randn(IdStruct.NYAMU,1)*AMUWNVar.^0.5;
    
            meas(IdStruct.IdYAMU) = meas(IdStruct.IdYAMU) + AMUWN + AMUBias;
        end
        if IdStruct.NYGyro
            % gyrownstddev = randn(IdStruct.NYSST,1)*SimParam.GyroWNVar.^0.5;
            % meas(IdStruct.IdYSST) = meas(IdStruct.IdYSST) + gyrownstddev;
            GyroWNVar   = ModelStruct.Sens.GyroWNVar;
            GyroBias    = ModelStruct.Sens.GyroBias;
            GyroWN      = randn(IdStruct.NYGyro,1)*GyroWNVar.^0.5;
    
            meas(IdStruct.IdYGyro) = meas(IdStruct.IdYGyro) + GyroWN + GyroBias;
        end
        if IdStruct.NYRW
            % SimParam.RWWNVar
            RWWNVar   = ModelStruct.Sens.RWWNVar;
            RWWN = randn(IdStruct.NYRW,1)*RWWNVar.^0.5;
            
            meas(IdStruct.IdYRW) = meas(IdStruct.IdYRW) + RWWN;
        end
        if IdStruct.NYSST
            % white noise
            SSTWNVar = ModelStruct.Sens.SSTWNVar;
            SSTWN    = randn(IdStruct.NYSST,1)*SSTWNVar.^0.5;
    
            meas(IdStruct.IdYSST) = meas(IdStruct.IdYSST) + SSTWN;
        end
    end
end

