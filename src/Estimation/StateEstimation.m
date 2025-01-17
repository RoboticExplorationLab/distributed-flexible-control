function [x_hat,y_hat] = StateEstimation(EstStruct, Meas, CtrlComm, wPert,...
                          nom_x, x_hat_prev, G, SimParam, istep)
    % Navigation/estimation function. State estimation approach is
    % different for the different states/measurements of the system.
    A  = EstStruct.A;
    B  = EstStruct.B;
    C  = EstStruct.C;
    D  = EstStruct.D;
    ikf = 1+mod(istep-1,SimParam.NStepSST);
    L  = squeeze(EstStruct.L(:,:,ikf));
    
    % linear state around operating point
    dx_hat_prev = x_hat_prev - nom_x;

    % y_hat_prev = C*x_hat_prev + D*CtrlComm; % luen-like
    dx_hat      = A*dx_hat_prev + B*CtrlComm + G*wPert;
    y_hat_prev  = C*(nom_x + dx_hat) + D*CtrlComm; % kf-like
    res = (Meas - y_hat_prev);
    dx_hat = dx_hat + L*res;

    % recovering state
    x_hat = nom_x + dx_hat;

    % measurement post-residual
    y_hat = C*x_hat + D*CtrlComm;
end

