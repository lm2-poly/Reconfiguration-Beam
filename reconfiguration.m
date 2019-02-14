%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM reconfiguration.m
% THE THIRD ORDER DIFFERENTIAL EQUATION CODED IN beamwithdrag.m IS
% INTEGRATED NUMERICALLY WITH A RUNGE-KUTTA ALGORITHM. THE SHOOTING METHOD
% IS USED TO CONVERT THE BOUNDARY VALUE PROBLEM INTO AN INITIAL VALUE
% PROBLEM. THE ANGLE AT THE FREE END OF THE BEAM IS INITIALLY GUESSED AND A
% MULLER ALGORITHM IS USED TO ITERATE TO THE CORRECT VALUE OF THE FREE-END
% ANGLE.
%
% PROGRAM REQUIRES beamwithdrag.m
%
% MORE INFORMATION CAN BE FOUND IN THE PAPER "DRAG REDUCTION OF FLEXIBLE
% PLATES BY RECONFIGURATION" BY F. GOSSELIN, E DE LANGRE AND B.A.
% MACHADO-ALMEIDA. SEE ALSO THE THESIS BE F. GOSSELIN "MECANISMES
% D'INTERACTIONS FLUIDE-STRUCTURE ENTRE ECOULEMENTS ET VEGETATION."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact
format long
clear resudrag 

%TOLERANCES
tol=0.000000001;
tolrel=0.000001;
options = odeset('RelTol',10^(-5),'AbsTol',tol);
maxiter=200;

%INITIALISATION OF VARIABLES
lastt1=0;   %VALUE OF THE FREE-END ANGLE AT THE LAST CYCD VALUE
lastt2=0;   %VALUE OF THE FREE-END ANGLE AT THE SECOND TO LAST CYCD VALUE
lastt3=0;   %VALUE OF THE FREE-END ANGLE AT THE THIRD TO LAST CYCD VALUE

beta=1      %W_1/W_1 TAPERING RATIO (A VALUE OF 1 IS FOR A STRAIGTH BEAM)
NR=100      %NUMBER OF CYCD POINTS
cymin=-2;   %SMALLEST VALUE OF CYCD = 10^cymin
cymax=5;    %LARGEST VALUE OF CYCD = 10^cymax

for nres=1:NR       %CYCD LOOP
    CYCD=10^(log10(10^cymin)+(nres-1)/(NR-1)*(log10(10^cymax)-log10(10^cymin)));
                    %EVALUATE THE NEXT VALUE OF CYCD FOR EQUIDISTANT POINTS
                    %ON THE LOGARITHMIC PLOT
    
    %INITIAL GUESSES OF THE FREE-END ANGLE. MULLER ALGORITHM REQUIRES THREE
    %GUESSES.
    if nres==1          %FOR THE FIRST POINT THE GUESSES ARE NOTHING BUT GUESSES
        g0=1;
        g1=0.001;
        g2=pi/2;
    elseif nres==2      %FOR THE SECOND VALUE OF CYCD, WE USE THE VALUE OF THE SOLUTION AT THE PREVIOUS POINT
        g2=lastt1*1.01;
        g1=lastt1*0.99;
        g0=lastt1;
    elseif nres==3      %FOR THE THIRD VALUE OF CYCD, WE USE THE PREVIOUSLY FOUND SOLUTIONS
        g2=lastt1*0.99;
        g1=lastt2;
        g0=lastt1;
    else                %FOR ALL HIGHER VALUES OF CYCD, WE USE THE THREE PREVIOUSLY FOUND SOLUTIONS
        g2=lastt3;
        g1=lastt2;
        g0=lastt1;
    end

    %WE CALL THE RUNGE-KUTTA ALGORITHM TO INTEGRATE OUR DIFFERENTIAL
    %EQUATION "beamwithdrag" OVER THE LENGTH OF THE BEAM FROM THE FREE END
    %TO THE CLAMPED END. THE ALGORITHM RETURNS THE VALUES OF THETA, THETA'
    %AND THETA'' (IN THE "TH" ARRAY) IN FUNCTION OF THE LAGRANGIAN COORDINATE S.
    %THE INTEGRATION IS PERFORMED USING THE FIRST GUESS OF THE END ANGLE
    %g1.
    [S,TH]=ode45(@(s,th) beamwithdrag(s,th,CYCD,beta),[1,0],[g1,0,0]',options);
    [nn garbage]=size(TH);
    f1=TH(nn,1)-pi/2;   %THE ANGLE AT THE CLAMPED END SHOULD BE 0. THE VALUE f1 IS THUS A MEASURE OF THE ERROR

    %THE INTEGRATION IS PERFORMED USING THE SECOND GUESS OF THE END ANGLE
    %g2.
    [S,TH]=ode45(@(s,th) beamwithdrag(s,th,CYCD,beta),[1,0],[g2,0,0]',options);
    [nn garbage]=size(TH);
    f2=TH(nn,1)-pi/2;
    
    iter=1;         %INITIALISATION OF THE ITERATION COUNTER
    enough=0;       %FLAG
    er=10000*tol;   %INITIALISATION OF THE ABSOLUTE ERROR
    relaer=10000;   %INITIALISATION OF THE RELATIVE ERROR
    while (iter<maxiter+1)  %ITERATIVE LOOP TO CONVERGE TO THE CORRECT FREE-END ANGLE USING A MULLER ALGORITHM

        %THE INTEGRATION IS PERFORMED USING THE GUESS OF THE END ANGLE g0.
        %THE VALUE OF g0 IS UPDATED AT EVERY LOOP
        [S,TH]=ode45(@(s,th) beamwithdrag(s,th,CYCD,beta),[1,0],[g0,0,0]',options);
        [nn garbage]=size(TH);
        f0=TH(nn,1)-pi/2;
        
        %PRINTING THE LAST GUESS AND THE LAST ERROR
        status=strcat('lastguess=',num2str(g0,16'),', angleat0=',num2str(f0,16))
    
        %MULLER ALGORITHM TO FIND THE NEXT GUESS
        a0=f0;
        aa=(f1-f0)/(g1-g0);
        a2=((f2-f1)/(g2-g1)-aa)/(g2-g0);
        a1=aa+(g0-g1)*a2;

        g2=g1;
        f2=f1;
        g1=g0;
        f1=f0;

        %FROM MULLER, WE OBTAIN TO VALUES OF THE END-ANGLE. WE KEEP THE ONE
        %CLOSEST TO THE PREVIOUS GUESS.
        pg(1)=g0+2*a0/(-a1-sqrt(a1^2-4*a0*a2));
        pg(2)=g0+2*a0/(-a1+sqrt(a1^2-4*a0*a2));
        [garbage,n1]=min(abs(g0-pg));
        ppg=pg(n1);

        if n1==1
            n2=1;
        else
            n2=2;
        end

        er=abs(g0-ppg);
        if ppg==0
            relaer=1;
        else
            relaer=abs(er/ppg);
        end

        %CRITERIA FOR CONVERGENCE
        if (er<tol&relaer<tolrel)%&(fi<10^(-10))|((ri==2)&(relaer<100*tolrel)&(fi<10^(-15)))%(er<tol&relaer<tolrel*10) | (er<tol*10&relaer<tolrel)
            answer(1,1)=g0;
            answer(1,2)=f0;
            answer(1,3)=iter;
            answer(1,4)=er;
            answer(1,5)=relaer
            break
        end
        g0=ppg;
        iter=iter+1;
    end
    
    %UPDATING THE VALUE OF THE LAST CONVERGED VALUES
    lastt3=lastt2;
    lastt2=lastt1;
    lastt1=g0;
    
    if iter>maxiter
        WARNINGMAXITERREACHED=1
        break
    end
    
    %ONCE A CONVERGED VALUE OF THETA IS FOUND, THE INTEGRATION IS PERFOMED
    %ONE LAST TIME WITH MORE PRECISION TO OBTAIN A PRECISE FUNCTION OF THE
    %DEFORMATION
    npint=-10000;%min(ceil(kap)*5,-100);
    [S,TH]=ode45(@(s,th) beamwithdrag(s,th,CYCD,beta),[1:1/npint:0],[answer(1,1),0,0]',options);
    [nn garbage]=size(TH);
    
    %INTEGRATION OVER THE LENGTH OF THE BEAM OF THE DRAG FORCE
    dsi=-(S(nn,1)-S(nn-1,1));  %BEAM ELEMENT LENGTH
    clear drag
    drag=0;
    for n=nn:-1:1
        theta=TH(n,1);
        pos(n,1)=S(n,1);
        sn=pos(n,1);
        if n==nn
            pos(n,2)=0;
            pos(n,3)=0;
            drag=drag+2*(1+(beta-1)*sn)/(beta+1)*sin(theta)^3*dsi;
        else
            pos(n,2)=pos(n+1,2)+dsi*sin(theta);
            pos(n,3)=pos(n+1,3)-dsi*cos(theta);
            if n==1
                drag=drag+2*(1+(beta-1)*sn)/(beta+1)*sin(theta)^3*dsi/2;
            else
                drag=drag+2*(1+(beta-1)*sn)/(beta+1)*sin(theta)^3*dsi;      %INTEGRATION OF THE DRAG
            end
        end
    end
    
    %SAVING THE RESULTS IN AN ARRAY
    resudrag(nres,1)=CYCD;
    resudrag(nres,2)=drag;
    
    %PLOTTING THE LAST FOUND POINT ON THE RECONFIGURATION CURVE
    figure(11)
    loglog(CYCD,drag,'ok')
    hold on
    grid on
    xlabel('C_Y C_D')
    ylabel('R')
    
    %PLOTTING THE SHAPE OF THE DEFORMED BEAM
    figure(2)
    hold on
    plot(-pos(:,2),pos(:,3),'-k')
    plot(pos(:,2),pos(:,3),'-k')
    %axis equal
    set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[13,6,13])
    set(gca,'XTickLabel',{})
    set(gca,'YTickLabel',{})
    box off
    ylim([-1 0])
    xlim([-1.08 1.08])
    axis off
    
end

%UNCOMMENT THE NEXT LINE TO SAVE THE RESULTS ON THE HARD DRIVE
%csvwrite('basicmodel-dragVScycd-beta=1.csv',resudrag);
