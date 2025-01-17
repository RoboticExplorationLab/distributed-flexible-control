function PostProMain(ResStruct, SimParam, IdStruct, ModelStruct)
close all;

tVec      = ResStruct.tVec;
UVec      = ResStruct.UVec;
CLStates = ResStruct.CLStates;
OLStates = ResStruct.OLStates;
x_ref    = ResStruct.x_ref; 
u_ref    = ResStruct.u_ref;
y_hat    = ResStruct.y_hat;
CLMeas   = ResStruct.CLMeas;
NVib     = SimParam.PostPro.NVib;
minrwspeed  = ModelStruct.Act.minrwspeed;
nomrwspeed  = ModelStruct.Act.nomrwspeed;
maxrwspeed  = ModelStruct.Act.maxrwspeed;
MaxTorqueRW = ModelStruct.Act.MaxTorqueRW;
tdiff       = 100;
file        = SimParam.PostPro.plotfile;

figID = 1;
% % Estimation 
if SimParam.PostPro.SensPlots
    % Plot Measurements
    figID = figID+1;
    PlotMeasurements(tVec, CLMeas, IdStruct, figID, 'Real', file);
    % Plot Measurement Estimates
    figID = figID+1;
    PlotMeasurements(tVec, y_hat, IdStruct, figID, 'Est.', file);
    % Plot Measurement Estimate Error
    figID = figID+1;
    PlotMeasurementError(tVec, y_hat, CLMeas, IdStruct, figID, file);
end

% Control Plots
if SimParam.PostPro.ActPlots
    % % % Control
    % % RW Torques
    figID = figID+1;
    PlotRWTorques(tVec,u_ref, UVec,MaxTorqueRW,IdStruct,figID, file);
    % 
    % % RW Speed
    figID = figID+1;
    PlotRWSpeeds(tVec, x_ref, CLStates, nomrwspeed, minrwspeed, maxrwspeed,...
                                                    IdStruct, figID, file);
    % % RW Speed and torques
    figID = figID+1;
    PlotRWTorquesNice(tVec, u_ref, UVec, x_ref, CLStates, MaxTorqueRW,...
               nomrwspeed, minrwspeed, maxrwspeed, IdStruct, figID, file);
end

% % Rigid Body attitude
if SimParam.PostPro.att
    figID = figID+1;
    PlotAttitude(tVec,x_ref,CLStates,OLStates,IdStruct, figID, file);
    % angular velocity
    figID = figID+1;
    PlotAngVel(tVec,x_ref,CLStates,OLStates,IdStruct, figID, file);
    % Attitude Error
    figID = figID+1;
    PlotAttitudeError(tVec,x_ref, CLStates, OLStates, IdStruct, tdiff, figID, file);
end

if SimParam.PostPro.vibmodes
    % % Vibration modes
    figID = figID+1;
    PlotVibrationModes(tVec,CLStates,OLStates,IdStruct, figID, file);
    % 
    % Vibration Modes Separately FullMode
    figID = figID+1;
    PlotVibModesSep(tVec, CLStates, OLStates, IdStruct, NVib, ...
                                                            figID, file);
    % 
    % RMS Deformation
    figID = figID+1;
    PlotRMSDefErr(tVec, CLStates, OLStates, ModelStruct, IdStruct, ...
                                                      tdiff, figID, file);
end


    % local functions
    function PlotRWTorques(tVec,u_ref,UVec,MaxTorqueRW,IdStruct,figID, file)
        Nctrlpoints = IdStruct.Nctrlpoints;
        NT = length(tVec);
        fig = figure(figID);
        ylab = {'X [Nm]','Y [Nm]','Z [Nm]'};
        % assuming all RW are facing x,y,z at every node in order
        for j= 1:min(3,max(IdStruct.RWperpoint))
            subplot(min(3,max(IdStruct.RWperpoint)),1,j)
            hold on;
            plot(tVec,MaxTorqueRW(1)*ones(1,NT),'k--');
            plot(tVec,-MaxTorqueRW(1)*ones(1,NT),'k--');
            for i=1:length(IdStruct.RWperpoint)
                if IdStruct.RWperpoint(i) >= j
                    stairs(tVec,u_ref(sum(IdStruct.RWperpoint(1:i-1))+j,:),'k')
                    stairs(tVec,UVec(:,sum(IdStruct.RWperpoint(1:i-1))+j),...
                                                        'r','LineWidth',0.5)
                end
            end
            grid on;
            ylabel(ylab{j})
            if j == 1
                title('RW torque Commands')
            end
        end
        xlabel('time [s]')
        savefig(fig,[file 'RWTorques.fig']);
        saveas(fig,[file 'RWTorques.png']);
    end

    function PlotRWTorquesNice(tVec,u_ref,UVec,x_ref,CLStates,MaxTorqueRW,...
                nomrwspeed, minrwspeed, maxrwspeed, IdStruct, figID, file)
        NT = length(tVec);
        fig = figure(figID);
        subplot(2,1,1)
        hold on;
        plot(tVec,MaxTorqueRW(1)*ones(1,NT),'k--');
        
        for j= 1:IdStruct.NU
            stairs(tVec,u_ref(j,:),'k')
            stairs(tVec,UVec(:,j),'r','LineWidth',0.5)
        end
        plot(tVec,-MaxTorqueRW(1)*ones(1,NT),'k--');
        grid on;
        ylabel('$\tau_{RW}$ [Nm]','Interpreter','latex','fontsize',12)
        legend('RW Torque limit','RW Torque','fontsize',10)
        ylim(1.1*[-MaxTorqueRW(1),MaxTorqueRW(1)])
        xlim([0,20])
        xlabel('time [s]','fontsize',10)
        subplot(2,1,2)
        hold on;
        plot(tVec,(maxrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
        for i=1:IdStruct.RealMdl.NXRW
            plot(tVec,(x_ref(IdStruct.RealMdl.IdXRW(i),:))*60/(2*pi),'k')
            plot(tVec,(CLStates(:,IdStruct.RealMdl.IdXRW(i))-nomrwspeed(i))*60/(2*pi),'r',...
                                                           'LineWidth',0.5)
        end
        plot(tVec,(minrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
        grid on;
        xlim([0,20])
        ylabel('$\omega_{RW}$ [rpm]','Interpreter','latex','fontsize',12)
        legend('RW Speed limit','RW Speed','fontsize',10)
        xlabel('time [s]','fontsize',10)
        savefig(fig,[file 'RWTorquesSpeed.fig']);
        saveas(fig,[file 'RWTorquesSpeed.png']);
    end
    
    function PlotRWSpeeds(tVec,x_ref,CLStates,nomrwspeed,minrwspeed,...
                                            maxrwspeed,IdStruct,figID,file)
        Nctrlpoints = IdStruct.Nctrlpoints;
        Nrwpp = min(3,max(IdStruct.RWperpoint));
        NT = length(tVec);
        Nq = IdStruct.RealMdl.Nmod;
        fig=figure(figID);
        ylbls = {'X [rpm]','Y [rpm]','Z [rpm]'};
        for j=1:Nrwpp
            subplot(Nrwpp,1,j)
            subplot(Nrwpp,1,j)
            hold on;
            plot(tVec,(maxrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
            plot(tVec,(minrwspeed-nomrwspeed)*60/(2*pi)*ones(1,NT),'k--');
            
            for i=1:length(IdStruct.RWperpoint)
                Nrwm1 = sum(IdStruct.RWperpoint(1:i-1));
                if IdStruct.RWperpoint(i) >= j
                    wrw = CLStates(:,IdStruct.RealMdl.IdXRW(Nrwm1+j))*60/(2*pi);
                    dwrw = wrw - nomrwspeed(Nrwm1+j)*60/(2*pi);
                    plot(tVec,x_ref(IdStruct.RealMdl.IdXRW(Nrwm1+j),:)*60/(2*pi),'k');
                    plot(tVec,dwrw,'r','LineWidth',0.5);
                end
            end
            grid on;
            ylabel(ylbls{j})
            if j ==1
                title('RW Speed wrt nominal Speed')
            end
        end
        xlabel('time [s]','fontsize',12)
        savefig(fig,[file 'RWSpeeds.fig']);
        saveas(fig,[file 'RWSpeeds.png']);
    end
    
    function PlotVibrationModes(tVec,CLStates,OLStates,IdStruct, figID, file)
        Nctrlpoints = IdStruct.Nctrlpoints;
        Nq = IdStruct.RealMdl.Nmod;
        fig=figure(figID);
        subplot(2,2,1)
        hold on;
        for i=1:Nq
            plot(tVec,OLStates(:,IdStruct.RealMdl.IdXmodq(i)))
        end
        grid on;
        ylabel('eta')
        title('Vibration modes OL')
        subplot(2,2,2)
        hold on;
        for i=1:Nq
            plot(tVec,OLStates(:,IdStruct.RealMdl.IdXmodd(i)))
        end
        grid on;
        ylabel('etadot [Hz]')
        subplot(2,2,3)
        hold on;
        for i=1:Nq
            plot(tVec,CLStates(:,IdStruct.RealMdl.IdXmodq(i)))
        end
        grid on;
        title('Vibration modes CL')
        subplot(2,2,4)
        hold on;
        for i=1:Nq
            plot(tVec,CLStates(:,IdStruct.RealMdl.IdXmodd(i)))
        end
        grid on;
        xlabel('time [s]')
    end
    
    function PlotVibModesSep(tVec, CLStates, OLStates, IdStruct, Nmodes, figID, file)
        Nq = IdStruct.RealMdl.Nmod;
        CLModPos =  CLStates(:,IdStruct.RealMdl.IdXmodq); 
        OLModPos =  OLStates(:,IdStruct.RealMdl.IdXmodq);
        CLModVel =  CLStates(:,IdStruct.RealMdl.IdXmodd);
        OLModVel =  OLStates(:,IdStruct.RealMdl.IdXmodd);
        fig = figure(figID);
        ylab = {'X [Nm]','Y [Nm]','Z [Nm]'};
        for j= 1:Nmodes
            subplot(Nmodes,2,2*(j-1)+1)
            hold on;
            plot(tVec,CLModPos(:,j),'Color',[0,0,0])
            plot(tVec,OLModPos(:,j),'--','Color',[0,0,0])
            grid on;
            ylabel(['$\eta_' num2str(j) '$'],'Interpreter','latex','fontsize',12)
            legend('Closed-loop','Open-loop','fontsize',10)
            xlabel('time [s]','fontsize',10)
            subplot(Nmodes,2,2*j)
            hold on;
            plot(tVec,CLModVel(:,j),'Color',[0,0,0])
            plot(tVec,OLModVel(:,j),'--','Color',[0,0,0])
            grid on;
            ylabel(['$\dot{\eta}_' num2str(j) '$'],'Interpreter','latex','fontsize',12)
            if j == 1
                title('$$\dot\eta$$ [Hz]','Interpreter','latex')
                legend('Closed-loop','Open-loop','fontsize',10)
            elseif j == Nmodes
                xlabel('time [s]','fontsize',12)
            end
        end
        savefig(fig, [file 'ModesOLvsCLposter2.fig']);
        saveas(fig, [file 'ModesOLvsCLposter2.png']);
    end

    function PlotRMSDefErr(tVec, CLStates, OLStates, ModelStruct, IdStruct, tdiff, figID, file)
        IdXmodq = IdStruct.RealMdl.IdXmodq;
        mass    = ModelStruct.SCMass;
        JCL     = vecnorm(CLStates(:,IdXmodq),2,2)/(mass^0.5);
        JOL     = vecnorm(OLStates(:,IdXmodq),2,2)/(mass^0.5);

        IdConv = tdiff/(tVec(2)-tVec(1)):length(JCL);
        errcl = sum(JCL(IdConv).^2/length(IdConv)).^0.5;
        errol = sum(JOL(IdConv).^2/length(IdConv)).^0.5;
        disp(['RMS Deformation, CL: ' num2str(1e6*errcl) ', OL: ' num2str(1e6*errol) ' microm'])

        fig=figure(figID);
        plot(tVec,1e6*JCL,'k','LineWidth',1.5)
        hold on
        plot(tVec,1e6*JOL,'k--','LineWidth',1.5)
        grid on;
        ylabel('RMS Def. [\mu m]','fontsize',12)
        legend('Closed Loop','Open Loop','fontsize',10)
        title('RMS Deformation [\mu m]')
        xlabel('time [s]','fontsize',10)
        savefig(fig,[file 'RMSDefOLvsCL.fig']);
        saveas(fig,[file 'RMSDefOLvsCL.png']);
    end
    
    function PlotMeasurements(tVec, y_hat, IdStruct, figID, type, file)
        NYVec = [IdStruct.NYGyro,IdStruct.NYAMU,IdStruct.NYRW];
        IdYVec = {IdStruct.IdYGyro,IdStruct.IdYAMU,IdStruct.IdYRW};
        Nsp = sum(NYVec~=0);
        ylbls = {'Gyro [rad/s]','AMU [m/s^2]','RW [rad/s]'};
        
        k = 1;
        fig=figure(figID);
        for i=1:Nsp
            k = find(NYVec(k:end)~=0,1,'first')+k-1;
            idYi = IdYVec{k};
            NYi  =  NYVec(k);
            subplot(Nsp,1,i)
            hold on;
            for j=1:NYi
                plot(tVec,y_hat(idYi(j),:),'k')
            end
            grid on;
            ylabel(ylbls(k))
            if i==1
                title([type ' measurements'])
            end
            k = k + 1;
        end
        xlabel('time [s]');
        savefig(fig,[file 'TrueMeas.fig']);
        saveas(fig,[file 'TrueMeas.png']);
    end
    
    function PlotMeasurementError(tVec, y_hat, CLMeas, IdStruct, figID, file)
        NYVec = [IdStruct.NYGyro,IdStruct.NYAMU];%,IdStruct.NYRW];
        IdYVec = {IdStruct.IdYGyro,IdStruct.IdYAMU};%,IdStruct.IdYRW};
        Nsp = sum(NYVec~=0);
        ylbls = {'Gyro Est. resid. [rad/s]','AMU [m/s^2]','RW [rad/s]'};
        k = 1;
        fig=figure(figID);
        for i=1:Nsp
            k = find(NYVec(k:end)~=0,1,'first')+k-1;
            idYi = IdYVec{k};
            NYi  =  NYVec(k);
            subplot(Nsp,1,i)
            for j=1:NYi
                semilogy(tVec,abs(y_hat(idYi(j),:)-CLMeas(idYi(j),:)),'k')
                hold on;
            end
            grid on;
            ylabel(ylbls(k),'fontsize',12),%'Interpreter','latex'
            ylim([0,20])
            k = k + 1;
        end
        xlabel('time [s]','fontsize',10);
        savefig(fig,[file 'MeasErr.fig']);
        saveas(fig,[file 'MeasErr.png']);
    end
    
    function PlotAttitudeError(tVec, x_ref, CLStates, OLStates, IdStruct, tdiff, figID, file)
        IdXrigrq = IdStruct.RealMdl.IdXrigrq;
        IdXrigrd = IdStruct.RealMdl.IdXrigrd;
        NT      = length(tVec);
        quatref = mrp2quat(x_ref(IdXrigrq,:)');
        quatref(~any(quatref~=0,2),1) = 1;
        zol = quatmultiply(quatconj(quaternion(mrp2quat(OLStates(:,IdXrigrq)))),...
                           quaternion(quatref)); 
        zcl = quatmultiply(quatconj(quaternion(mrp2quat(CLStates(:,IdXrigrq)))),...
                           quaternion(quatref)); 
        zol = compact(zol);
        zcl = compact(zcl);

        OLAttErr = abs(2*atan2d(vecnorm(zol(:,2:4),2,2),zol(:,1))*3600);
        CLAttErr = abs(2*atan2d(vecnorm(zcl(:,2:4),2,2),zcl(:,1))*3600);
        % angle = angle_between_vecs(v1, v2);
        OLAngVelNorm = vecnorm(OLStates(:,IdXrigrd),2,2);
        CLAngVelNorm = vecnorm(CLStates(:,IdXrigrd)-x_ref(IdXrigrd,:)',2,2);
        % convergence time: once angvel reaches steady state (mean ang vel
        % within tol of its covar for a 1s int time)
        
        IdConv = tdiff/(tVec(2)-tVec(1)):length(CLAttErr);
        aperms = sum(CLAttErr(IdConv).^2/length(IdConv)).^0.5;
        disp(['Pointing APE RMS: ' num2str(aperms) ' arcsec'])
        
        fig=figure(figID);
        subplot(2,1,1)
        hold on;
        plot(tVec, OLAttErr,'r--')
        plot(tVec, CLAttErr,'r')
        grid on;
        xregion(0,tdiff,FaceColor='k',FaceAlpha=0.1,EdgeColor='k');
        xregion(tdiff,tVec(end),FaceColor='y',FaceAlpha=0.1,EdgeColor='k');
        plot([tdiff,tVec(end)],aperms*[1,1], '--k')
        ylabel('Point. Error [arcsec]','Interpreter','latex','fontsize',12)
        xlabel('time [s]','fontsize',10)
        legend('Open-loop','Closed-loop','fontsize',10)
        ylim([0, 1.5*max(CLAttErr)])
        xlim([tVec(1),tVec(end)])
        % title('$$Attitude Quaternion$$','Interpreter','latex')
        subplot(2,1,2)
        hold on;
        plot(tVec, rad2deg(OLAngVelNorm),'r--')
        plot(tVec, rad2deg(CLAngVelNorm),'r')
        grid on;
        ylabel('$$\|\omega\|$$ [deg/s]','Interpreter','latex','fontsize',12)
        xlabel('time [s]','fontsize',10)
        legend('Open-loop','Closed-loop','fontsize',10)
        ylim([0, 1.5*max(rad2deg([OLAngVelNorm, CLAngVelNorm]),[],'all')])
        xlim([tVec(1),tVec(end)])
        savefig(fig,[file 'AttitudeErr.fig']);
        saveas(fig,[file 'AttitudeErr.png']);
    end

   
    function PlotAttitude(tVec, x_ref, CLStates, OLStates, IdStruct, figID, file)
        % NXr = IdStruct.NXr;
        IdXrigrq = IdStruct.RealMdl.IdXrigrq;
        Nrigrq   = IdStruct.RealMdl.Nrigrq;
        fig=figure(figID);
        subplot(Nrigrq,1,1)
        hold on;
        plot(tVec, x_ref(IdXrigrq(1),:),'k')
        plot(tVec, OLStates(:,IdXrigrq(1)),'r--')
        plot(tVec, CLStates(:,IdXrigrq(1)),'r')
        grid on;
        ylabel('$$p_1$$','Interpreter','latex')
        title('$$Attitude Quaternion$$','Interpreter','latex')
        subplot(Nrigrq,1,2)
        hold on;
        plot(tVec, x_ref(IdXrigrq(2),:),'k')
        plot(tVec, OLStates(:,IdXrigrq(2)),'r--')
        plot(tVec, CLStates(:,IdXrigrq(2)),'r')
        grid on;
        ylabel('$$p_2$$','Interpreter','latex')
        subplot(Nrigrq,1,3)
        hold on;
        plot(tVec, x_ref(IdXrigrq(3),:),'k')
        plot(tVec, OLStates(:,IdXrigrq(3)),'r--')
        plot(tVec, CLStates(:,IdXrigrq(3)),'r')
        grid on;
        ylabel('$$p_3$$','Interpreter','latex')
        xlabel('time [s]')
        savefig(fig,[file 'Attitude.fig']);
        saveas(fig,[file 'Attitude.png']);
    end

    function PlotAngVel(tVec, x_ref, CLStates, OLStates, IdStruct, figID, file)
        % NXr = IdStruct.NXr;
        IdXrigrd = IdStruct.RealMdl.IdXrigrd;
        Nrigrd   = IdStruct.RealMdl.Nrigrd;
        fig=figure(figID);
        subplot(Nrigrd,1,1)
        hold on;
        plot(tVec, x_ref(IdXrigrd(1),:),'k')
        plot(tVec, CLStates(:,IdXrigrd(1)),'r')
        plot(tVec, OLStates(:,IdXrigrd(1)),'r--')
        grid on;
        ylabel('$$\omega_x$$ [rad/s]','Interpreter','latex')
        title('Angular velocity','Interpreter','latex')
        subplot(Nrigrd,1,2)
        hold on;
        plot(tVec, x_ref(IdXrigrd(2),:),'k')
        plot(tVec, OLStates(:,IdXrigrd(2)),'r--')
        plot(tVec, CLStates(:,IdXrigrd(2)),'r')
        grid on;
        ylabel('$$\omega_y$$ [rad/s]','Interpreter','latex')
        subplot(Nrigrd,1,3)
        hold on;
        plot(tVec, x_ref(IdXrigrd(3),:),'k')
        plot(tVec, OLStates(:,IdXrigrd(3)),'r--')
        plot(tVec, CLStates(:,IdXrigrd(3)),'r')
        grid on;
        ylabel('$$\omega_z$$ [rad/s]','Interpreter','latex')
        xlabel('time [s]')
        savefig(fig,[file 'AngVel.fig']);
        saveas(fig,[file 'AngVel.png']);
    end
end
