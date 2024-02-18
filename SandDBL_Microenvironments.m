
velocity = 0.1;

bottom_o2 = 200;
bottom_NO3 = 5;

grain_size=700;
kinVis = 1.238E-06;
g = 9.81;
por = 0.40;
permeability = 9.869*10^-13*735*10^6*(grain_size*1e-6)^2;

wavelength = 490*grain_size*1e-6;   % Relationship Wavelength

nitrification = 0;
anaerobicDen = 1;
o2_consumption = 40; %anaerobicDen*10;
aerobicDen = 0; %anaerobicDen;

conductivity = permeability ./ kinVis .* g;
k = 2*pi/wavelength;
dynamicHead = 1000*velocity^2*0.1/10000;
porewaterVel = dynamicHead*k*conductivity;
%porewaterVel = 10e-6;

q = 24*3600*k*conductivity*dynamicHead/pi/wavelength;  % Volume flux in m3 / m2 / s

timeO2 = bottom_o2./o2_consumption;
m_o2 = 1/k*log(0.42*k^2*conductivity*dynamicHead*timeO2*3600/por+1); % oxygen penetration depth
Flux_o2 = m_o2 * o2_consumption*24;
singleCell = 955;%*(o2_consumption/140);

%Anoxia Calculation
layers=4;
Vol=zeros(layers,1);
for n = 1:layers
    if n == layers
        m = 3.9;
    else
        m = n;
    end
    
    time = timeO2./layers*m;

    o2Conc = bottom_o2-time*o2_consumption;
    sandDBL = (grain_size*1e-6/(1+(0.62*((porewaterVel*grain_size*1e-6/2/kinVis)^0.41)*(1000 ^0.33))))^2/0.000000001*singleCell*1000/3600/o2Conc;
    sandDBL
    if sandDBL < 10
        Vol(n) = 0;
    else
        Vol(n) = 1 -(1.87* (sandDBL^(-0.26)));
    end
    if sandDBL > 1000
        Vol(n) = 1;
    end

end

anoxia = nanmean(Vol);

interface_NO3 = bottom_NO3-(anoxia*anaerobicDen)*timeO2;
timeNO3 = (bottom_NO3+nitrification*timeO2)./anaerobicDen;
m_NO3 = 1/k*log(0.42*k^2*conductivity*dynamicHead*(timeNO3+timeO2)*3600/por+1);
Flux_NO3 = ((aerobicDen*m_o2)*(1-anoxia)+(anaerobicDen*(m_NO3-m_o2)))*24
Flux_NO3_microEnv = m_o2*anoxia*anaerobicDen*24*0.4
Flux_NO3_microEnv./(Flux_NO3_microEnv+Flux_NO3)
