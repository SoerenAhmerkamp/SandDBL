clear all

%Fixed Parameters
bottom_o2 = 200;
bottom_NO3 = 6;
kinVis = 1.238E-06;
g = 9.81;
por = 0.37;

velocity = 0.1;
o2_consumption = 149;

[o2_consumptionGr, grain_sizeGr] = meshgrid(linspace(5,200,100),linspace(100,700,100));
ratio = zeros(size(o2_consumptionGr));
ano = zeros(size(o2_consumptionGr));

p = 0;
for xP = 1:size(o2_consumptionGr,1)
    for yO = 1:size(o2_consumptionGr,2)
        
    o2_consumption = o2_consumptionGr(xP,yO);
    grain_size = grain_sizeGr(xP,yO);
    
    %velocity = 0.1;
    
    p = p+1;
    permeability = 9.869e-13*735*(grain_size)^2*1e-6;
    %dynamicHead = 1.19e-3;
    wavelength = 490*grain_size*1e-6;

    anaerobicDen = o2_consumption*0.1;
    %nitrification = 5.8;
    %anaerobicDen = 2.57;
    %aerobicDen = 1.31;
    
    dynamicHead = 1000*velocity^2*0.1/10000; % 0.1 is the non-dimensional pressure head
    conductivity = permeability ./ kinVis .* g;
    k = 2*pi/wavelength;
    porewaterVel = dynamicHead*k*conductivity;

    timeO2 = bottom_o2./o2_consumption; % characteristic timescale O2
    timeNO3 = (bottom_NO3)./anaerobicDen; % charactersitic timescale NOx

    m_o2 = 1/k*log(0.42*k^2*conductivity*dynamicHead*timeO2*3600/por+1); 
    m_NO3 = 1/k*log(0.42*k^2*conductivity*dynamicHead*(timeNO3+timeO2)*3600/por+1);

    Flux_o2 = m_o2 * o2_consumption*24;
    singleCell = 955.*o2_consumption;
    %Anoxia Calculation
    layers=4;
    Vol=zeros(layers,1);
    for n = 1:layers
        if n == layers
            m = layers-0.1;
        else
            m = n;
        end
        
        time = timeO2./layers*m;
        m = 1/k*log(0.42*k^2*conductivity*dynamicHead*time*3600/por+1);
        porewaterVelD = porewaterVel * exp(-k*m);
        
        o2Conc = bottom_o2-time*o2_consumption;
        sandDBL = (grain_size*0.000001/(1+(0.62*((porewaterVelD*grain_size*0.000001/2/kinVis)^0.41)*(1000 ^0.33))))^2/0.000000001*singleCell*1000/3600/o2Conc;
        Vol(n) = 1 -(1.87* (sandDBL^(-0.26)));
    end

    anoxia = nanmean(Vol);
    Flux_NO3 = ((anaerobicDen*(m_NO3-m_o2)))*24;
    if anoxia < 0 | sandDBL < 30
       anoxia = 0; 
    end
    Flux_NO3_microEnv = m_o2*anoxia*anaerobicDen*24;
        
    ratio(xP,yO) = Flux_NO3_microEnv ./ (Flux_NO3+Flux_NO3_microEnv);
    ano(xP,yO) = anoxia;
    
    end
end

subplot(1,2,1)
contourf(grain_sizeGr,o2_consumptionGr, ano,30), shading flat, colorbar, axis square
ylabel('Consumption rate (µmol l^{-1} h^{-1})')
xlabel('Grain Size (µm)')
caxis([0, 1])
title('Anoxic Volumes')

subplot(1,2,2)
contourf(grain_sizeGr,o2_consumptionGr, ratio,30), shading flat, colorbar, axis square, colormap(flipud(hot()))
ylabel('Consumption rate (µmol l^{-1} h^{-1})')
xlabel('Grain Size (µm)')
caxis([0, 1])
title('Ratio oxic den / tot denitrification')

