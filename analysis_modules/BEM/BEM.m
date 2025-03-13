% Calculation of propeller characteristics based on the Blade Element
% Method model of L.L.M. Veldhuis
%
% Version 1.0 
% E.C.R. van Berkel (2012)
% Delft University of Technology
%
% -------------------------COMMENTS EvB----------------------------------
% Het programma is te runnen door ‘BEM.m’ te openen en te runnen.
% De propeller coefficienten die berekend worden, worden opgeslagen
% in de map ‘coefficienten’  met de bestandsnaam die verwijst naar V, J
% en ap (propeller invalshoek) bijvoorbeeld ‘coefficientsV30a0J0.75.mat’.
% 
% Er kunnen nu sweeps gedaan worden met snelheid, advance ratio, invalshoek
% van de propeller, en alle combinaties hiervan en dat is vrij gemakkelijk
% uit te breiden als andere mensen iets anders willen. Als er sweeps gedaan
% worden is het te adviseren om de plots en prints (regel 8 t/m 13) uit te
% zetten door de 1tjes te veranderen in 0. Als deze allemaal uit
% staan duurt het uitrekenen van de coefficienten ruim 1 seconde per
% datapunt.
% 
% De resultaten zijn hetzelfde als die van het Mathcad model.
% Ook dezelfde limitaties zitten er in. Het Mathcad model crashte
% bij hoge J en ap (bijvoorbeeld J = 1.5, ap = 8) en dat doet deze ook.
% Er staat namelijk een wortel in regel 138 en de inhoud van de
% wortel wordt dan negatief. Vervolgens krijgt onder andere de invalshoek
% een complexe waarde en doen de interpolaties van Cl en Cd het niet meer.
%-------------------------------------------------------------------------

clc;clear all;clf;close all;tic;

% State what needs to be displayed
plotBladeGeometry  = 0; % plot blade angle and blade chord
plotPolars         = 0; % plot the airfoils Cl and Cd vs aoa
plotLoopOutput     = 1; % plot phi, total velocity, Reynolds number and lift and drag, as well as axial and tangential inflow factor along the blade
plotFinalVariables = 0; % plot axial and swirl velocities, pressure rise, swirl angle, distributions of Cl and Cd along the blade and Normal force vs psi
printDistributions = 1; % print a, cl, cd, alast_ax and alast_tang distributions along the blade
printCoefficients  = 1; % print propeller coefficients (Tc, CT, Qc, Pc, Cp, CNp and n)

% Inputfiles
propdata          = importdata('propellerdata2.dat',' '); % contains propeller geometry
airfoildata1      = importdata('airfoil.dat',' ');       % contains airfoils polars

% Input of propeller parameters
B       = 6;      % number of blades
Dp      = .3048;  % propeller diameter            [m]
Beta075 = 36;     % blade angle at 3/4 R postion  [deg]
ducted  = 0;      % set to 1 if ducted, 0 if unducted

%Scale for drag coefficient to see its effect (1.0 means drag is unaltered)
factor_CD=1.00;

for U0  = 30      % Airspeed                      [m/s]
for J   = 1.0    % Advance Ratio
for ap  = 0       % effective propeller anlge of attack with respect to the incoming flow [deg]

% Constants
TC      = 15;    % air temperature               [deg celsius]
rho     = 1.225; % air density                   [kg/m3]

%% read the proller geometry file
nsections = length(propdata(:,1)); % number of input blade sections  
rR        = propdata(:,1);  % r/R, radial coordinate
cR        = propdata(:,2);  % c/R, chord normalized with blade length
beta      = propdata(:,3);  % beta with respect to the 3/4R position [deg]

if plotBladeGeometry == 1
figure()
subplot(2,1,1); plotEdited(rR,beta,'r',2,'r/R','beta (deg)','Blade angle distribution','on');ylim([1.1*min(beta) 1.1*max(beta)]);
subplot(2,1,2); plotEdited(rR,cR,'r',2,'r/R','c/R','Blade chord distribution','on');ylim([1.1*min(cR) 1.1*max(cR)]);
end

betarad     = (beta+Beta075)*pi/180;    % beta [rad]
aprad       = ap*pi/180;                % aoa propeller [rad]

%% Initial calculations
[Rp Sp q rps rpm r omega TK TK0 mu0 mu] = initcalc(rR,Dp,rho,U0,J,TC);

%% Determine number of stations to calculate
psinumber = 12;                      % number of stations to calculate
j         = 1:13;                    % psi counter
psi       = (360/psinumber) * (j-1); % psi vector [deg]
psirad    = psi*pi/180;              % psi vector [rad]

%% read and interpolate airfoil data
alfa = airfoildata1(:,1);
cl   = airfoildata1(:,2);
cdr   = airfoildata1(:,3);

%scale the drag coefficient now:
cdr=cdr*factor_CD;

x = -6:2:12;

clfit = fit(alfa,cl,'splineinterp');
cdfit = fit(alfa,cdr,'splineinterp');

if plotPolars == 1
figure()
subplot(2,1,1); plotEdited(x,feval(clfit,x),'b',2,'alpha (deg)','cl','Airfoil lift coefficient','on');hold on;
plot(alfa,cl,'o','LineWidth',2);legend('fit','data');legend('Location',[0.75 .65 0.1 .1]);
subplot(2,1,2); plotEdited(x,feval(cdfit,x),'b',2,'alpha (deg)','cd','Airfoil drag coefficient','on');hold on;
plot(alfa,cdr,'o','LineWidth',2);legend('fit','data');legend('Location',[0.75 .3 0.1 .1]);ylim([0 .2]);
end

%% Iteration loop
% initiate velocities for loop
va  = zeros(nsections,length(j)); % axial velocity                   [m/s]
vt  = zeros(nsections,length(j)); % tangential velocity              [m/s]

% induced velocities (set to zero for now)
vai = zeros(nsections,length(j)); % axial induced velocity           [m/s]
vti = zeros(nsections,length(j)); % tangential induced  velocity     [m/s]

vanew = zeros(nsections,length(j));

%% perform velocities loop
count = 0;
vdiff = 1;
pgcount = 0;
while vdiff > 0.005
    count = count + 1;
    for i = 1:nsections % i = spanwise section
        for j = 1:length(psirad) % j = psi, the blade angle in propeller plane 
            % calculate the velocity in the propeller plane
            vplane(i,j) = omega*r(i) - vt(i,j) + U0*sin(aprad)*sin(psirad(j));
            vplane(i,j) = vplane(i,j) + vti(i,j); % add tang. ind. vel.

            % calculate the angle phi for all blade sections: 
            phi(i,j)    = atan( ((U0+va(i,j))*cos(aprad)+vai(i,j)) / vplane(i,j));
            
            % calculate the local blade aoa
            a(i,j)      = betarad(i) - phi(i,j); 

            % Total velocity at the blade section
            w(i,j)      = sqrt( (U0*cos(aprad) + vai(i,j) + va(i,j))^2 + vplane(i,j)^2);

            % Mach number at the blade section
            Mach(i,j)   = w(i,j) / sqrt(1.4*287.05*TK);

            % Reynolds number
            Re(i,j)     = (rho*w(i,j)*cR(i)*Rp) / mu;

            % Determine Cl and Cd for Blade section
            Cl(i,j)     = feval(clfit,a(i,j)*180/pi);
            Cdr(i,j)     = feval(cdfit,a(i,j)*180/pi);

            % Apply Prandtl-Glauert correction when needed
            if Mach(i,j) > 0.7
                PG(i,j) = 1 / sqrt(1-(Mach(i,j)^2));
                Cl(i,j) = Cl(i,j) * PG(i,j); % apply PG correction5
                Cdr(i,j) = Cdr(i,j) * PG(i,j);
                pgcount = pgcount + 1;
            end

            % Calculate the loading
            dtdr(i,j)   = .5*rho*w(i,j)^2*B*cR(i)*Rp*(Cl(i,j)*cos(phi(i,j))-Cdr(i,j)*sin(phi(i,j))); % thrust component
            dqdr(i,j)   = .5*rho*w(i,j)^2*B*cR(i)*Rp*r(i)*(Cl(i,j)*sin(phi(i,j))+Cdr(i,j)*cos(phi(i,j))); % torque component
            dndr(i,j)   = dqdr(i,j)*sin(psirad(j)) / r(i); % normal force component

            phi1(i,j)   = asin( ((U0+va(i,j))*cos(aprad)+vai(i,j)) / w(i,j) ); % don't know what it is and what its use is
            phit(i,j)   = atan(rR(i)*tan(phi(i,j))); % shouldn't this be phi1? Seems to be no difference in outcome
            ss(i,j)     = (pi*Dp*sin(phit(i,j))) / B;
            
            % Prandtl Tip Loss factor
            if ducted ~= 1
                factorf(i,j)= (2/pi)*acos(exp(-pi*((Rp-r(i))/ss(i,j)))); 
            else
                factorf(i,j)= 1;
            end
% Next statement were changed to prevent sqrt of negative value. Be
% carefull when interpreting the results
            term_in_sqrt= (U0*cos(aprad)+vai(i,j))^2 + (dtdr(i,j) / (factorf(i,j)*pi*rho*r(i)))
            vanew(i,j)  = -.5*(U0*cos(aprad)+vai(i,j)) + .5*sqrt(term_in_sqrt);
            vtnew(i,j)  = dqdr(i,j) / (4*(pi*factorf(i,j)*rho*(r(i)^2)*(U0*cos(aprad)+vti(i,j)+vanew(i,j))));

            % Use relaxation
            va(i,j)     = .5*va(i,j) + .5*vanew(i,j);
            vt(i,j)     = .5*vt(i,j) + .5*vtnew(i,j);
            
            if count == 1
            a1ax        = va / (U0*cos(aprad));
            a1tang(i,j) = vt(i,j) / (omega*r(i));
            end
            
            alasttang(i,j) = vt(i,j) / (omega*r(i));
        end
    end
    vtdiff = max([max(max(vt-vtnew)) abs(min(min(vt-vtnew)))]);
    vadiff = max([max(max(va-vanew)) abs(min(min(va-vanew)))]);
    vdiff  = max([vtdiff vadiff]);
end
alastax = va / (U0*cos(aprad));

%% Calculate Propeller coefficients

for j = 1:length(psirad) % j = psi, the blade angle in propeller plane     
    for i = 1:nsections % i = spanwise section
        % swirl angle (deg)
        dswi(i,j) = (180/pi) * atan( (2*vt(i,j)*factorf(i,j)-U0*sin(aprad)*sin(psirad(j))) / (U0*cos(aprad)+factorf(i,j)*va(i,j)));

        % static pressure rise
        dps(i,j) = 2*rho*((omega*r(i))-(vt(i,j)*factorf(i,j)))*vt(i,j)*factorf(i,j);

        % total pressure rise
        dpt(i,j) = dps(i,j)+.5*rho*(2*vt(i,j)*factorf(i,j))^2;

        % dimenstionless total pressure rise
        dptq(i,j) = dpt(i,j) / q;
    end
    % Thrust
    tblade(j) = sum(((dtdr(2:end,j)+dtdr(1:end-1,j))/2).*(r(2:end)-r(1:end-1)));
    tblade(j) = tblade(j) + ((dtdr(end,j)/2)*(Rp-r(end)));
    
    % Torque
    qblade(j) = sum(((dqdr(2:end,j)+dqdr(1:end-1,j))/2).*(r(2:end)-r(1:end-1)));
    qblade(j) = qblade(j) + ((dqdr(end,j)/2)*(Rp-r(end)));

    % Normal Force
    nblade(j) = sum(((dndr(2:end,j)+dndr(1:end-1,j))/2).*(r(2:end)-r(1:end-1)));
    nblade(j) = nblade(j) + ((dndr(end,j)/2)*(Rp-r(end)));  
end
% Total values integrated over propeller
T = (1/(2*pi))* sum(((tblade(2:end)+tblade(1:end-1))/2).*(psirad(2:end)-psirad(1:end-1)));
Q = (1/(2*pi))* sum(((qblade(2:end)+qblade(1:end-1))/2).*(psirad(2:end)-psirad(1:end-1)));
FN = (1/(2*pi))* sum(((nblade(2:end)+nblade(1:end-1))/2).*(psirad(2:end)-psirad(1:end-1)));

% Coefficients
Tc  = T / (2*q*Dp^2);
CT  = T / (rho*rps^2*Dp^4);
Qc  = Q / (2*q*Dp^3);
Pc  = 2*pi*Qc/J;
Cp  = (2*pi*rps*Q) / (rho*rps^3*Dp^5);
CNp = FN / (q*Sp);
n   = Tc / Pc;


%% Display results
if pgcount > 0
    disp('Warning: Prandtl-Glauert corrections applied');
end

dirinfo = what;                     % provides info on current directory
curdir  = dirinfo.path;             % gives pathname of current directory
cd(strcat(curdir,'/coefficients')); % change directory to coefficients folder

%filename = num2str([U0 ap J],'coefficientsV%1.0fa%1.0fJ%1.2f.mat');
filename = num2str([U0 ap J],'coefficientsV%1.0fa%1.0fJ%1.2f.dat');

save(filename, 'Tc', 'CT', 'Qc', 'Pc', 'Cp', 'CNp', 'n'); % save propeller coefficients to MAT-file
cd ..;                              % change directory back to parent folder

if plotLoopOutput == 1
    figure()
    subplot(2,2,1); plotEdited(rR,phi(:,1)*180/pi,'r',2,'r/R','phi (deg)','Phi angle distribution','on');
    subplot(2,2,2); plotEdited(rR,w(:,1),'r',2,'r/R','w','Total velocity distribution','on');
    subplot(2,2,3); plotEdited(rR,Re(:,1),'r',2,'r/R','Re','Reynolds number ','on');
    subplot(2,2,4); plotEdited(rR,Cl(:,1),'r',2,'r/R','','Lift and drag polar','on');hold on;
    subplot(2,2,4); plot(rR,Cdr(:,1),'b','LineWidth',2);hold off;legend('Cl','Cd');legend('Location',[0.75 .15 0.1 .1]);

    figure()
    subplot(2,1,1); plotEdited(rR,a1ax(:,1),'r',2,'r/R','','Axial inflow factor','on');hold on;
    subplot(2,1,1); plotEdited(rR,alastax(:,1),'b',2,'r/R','','Axial inflow factor','on');legend('a1_a_x','alast_a_x');hold off;legend('Location',[0.75 .65 0.1 .1]);
    subplot(2,1,2); plotEdited(rR,a1tang(:,1),'r',2,'r/R','','Tangential inflow factor','on');hold on;
    subplot(2,1,2); plotEdited(rR,alasttang(:,1),'b',2,'r/R','','Tangential inflow factor','on');legend('a1_t_a_n_g','alast_t_a_n_g');hold off;legend('Location',[0.75 .15 0.1 .1]);
end

if printDistributions == 1
    disp('    rR        a(deg)    cl        cd        alast_ax  alast_tang');   
    disp('----------------------------------------------------------------');
    disp([rR a(:,1)*180/pi feval(clfit,a(:,1)*180/pi) feval(cdfit,a(:,1)*180/pi) alastax(:,1) alasttang(:,1)]);
end

if printCoefficients  == 1
    disp('');
    disp(num2str([U0 ap J],'V = %1.0f m/s, ap = %1.0f, J = %1.2f:'));
    disp('    Tc        CT        Qc        Pc        Cp        CNp       n');   
    disp('----------------------------------------------------------------------');
    disp([Tc CT Qc Pc Cp CNp n]);
    disp('');
end

if plotFinalVariables == 1   
    figure()
    subplot(2,2,1); plotEdited(rR,vt(:,4)/U0,'b',2,'','','','on');hold on;
    subplot(2,2,1); plotEdited(rR,va(:,4)/U0,'r',2,'','','','on');hold on;
    subplot(2,2,1); plotEdited(-rR,vt(:,10)/U0,'b',2,'','','','on');hold on;
    subplot(2,2,1); plotEdited(-rR,va(:,10)/U0,'r',2,'r/R','induced velocity','Axial and swirl velocities','on');hold off;
  
    subplot(2,2,2); plotEdited(rR,dps(:,4),'b',2,'','','','on');hold on;
    subplot(2,2,2); plotEdited(rR,dpt(:,4),'r',2,'','','','on');hold on;
    subplot(2,2,2); plotEdited(-rR,dps(:,10),'b',2,'','','','on');hold on;
    subplot(2,2,2); plotEdited(-rR,dpt(:,10),'r',2,'r/R','pressure rise','Total and static pressure rise (Pa)','on');hold off;
    
    subplot(2,2,3); plotEdited(rR,dpt(:,4)/q,'r',2,'','','','on');hold on;
    subplot(2,2,3); plotEdited(-rR,dpt(:,10)/q,'r',2,'r/R','pressure rise','Dimenstionless total pressure rise','on');hold off;
    
    subplot(2,2,4); plotEdited(rR,dswi(:,4),'r',2,'','','','on');hold on;
    subplot(2,2,4); plotEdited(-rR,dswi(:,10),'r',2,'r/R','swirl angle (deg)','Swirl angle (deg)','on');hold off;
    
    figure()
    subplot(3,1,1); plotEdited(rR,Cl(:,4),'r',2,'','','','on');hold on;
    subplot(3,1,1); plotEdited(-rR,Cl(:,10),'r',2,'r/R','Cl','Distribution of lift coefficient','on');hold off;
    subplot(3,1,2); plotEdited(rR,Cdr(:,4),'r',2,'','','','on');hold on;
    subplot(3,1,2); plotEdited(-rR,Cdr(:,10),'r',2,'r/R','Cd','Distribution of drag coefficient','on');hold off;
    subplot(3,1,3); plotEdited(psi,nblade,'r',2,'psi (deg)','nblade','Blade normal force vs. psi','on');xlim([0 360]);
end

end
end
end
disp('Calculations finished');
toc