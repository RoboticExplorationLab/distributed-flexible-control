function [xkp1,thetarw] = rigflexdynDisc(xk,uk, thetarw, RealStruct, IdStruct, ...
                                                    ModelStruct, SimParam)
    RealStruct.FiltFlag = true;
    DynFunc = @(t,x) rigflexcontDisc(x,uk, RealStruct, ...
                                        IdStruct, ModelStruct, SimParam);
    % RK4
    h   = RealStruct.SampTime; 
    xk  = [xk;thetarw];
    k_1 = DynFunc(0,xk);
    k_2 = DynFunc(0.5*h,xk+0.5*h*k_1); % DynFunc(xk+k_3*h)
    k_3 = DynFunc(0.5*h,xk+0.5*h*k_2);
    k_4 = DynFunc(h,xk+k_3*h);
    xkp1    = xk + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    thetarw = xkp1(IdStruct.RealMdl.NX+1:end);
    xkp1    = xkp1(1:IdStruct.RealMdl.NX);
end