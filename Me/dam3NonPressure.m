Nx = 32;
Ny = 32;
dx = 1 / Nx;
dt = 0.01;
gridDensity = zeros(Nx,Ny);
gridVelocityU = zeros(Nx,Ny);
gridVelocityV = zeros(Nx,Ny);
gridMass = zeros(Nx,Ny);

particleCount = 60;
particleVelocityU = zeros(particleCount,1);
particleVelocityV = zeros(particleCount,1);
particlePosX = zeros(particleCount,1);
particlePosY = zeros(particleCount,1);
particleCx = zeros(particleCount,1);
particleCy = zeros(particleCount,1);

for i = 1:particleCount
    particlePosX(i) = floor((i-1) / 10) * dx;
    particlePosY(i) = mod(i-1,10) * dx;
    particleVelocityU(i) = (0.5 - particlePosY(i));
    particleVelocityV(i) = -1;
end
time = 0;
timeFinal = 10;
while(time < timeFinal)
    gridVelocityU = zeros(Nx,Ny);
    gridVelocityV = zeros(Nx,Ny);
    gridWeight = zeros(Nx,Ny);
    % particle To Grid Velocity
    for k = 1:particleCount
        
        x = max(particlePosX(k) / dx + 1,1);
        y = max(particlePosY(k) / dx + 1,1);
        gridx = floor(x);
        gridy = floor(y);
        fx = x - gridx;
        fy = y - gridy;
        wx = [0.5 * (1.5 - fx).^2 0.75 - (fx - 1).^2 0.5*(fx - 0.5).^2];
        wy = [0.5 * (1.5 - fy).^2 0.75 - (fy - 1).^2 0.5*(fy - 0.5).^2];
        for j = 0:2
            for i = 0:2
                weight = wx(i+1) * wy(j+1);
                gridVelocityU(gridx+i,gridy+j) = gridVelocityU(gridx+i,gridy+j) + weight * particleVelocityU(k);
                gridVelocityV(gridx+i,gridy+j) = gridVelocityV(gridx+i,gridy+j) + weight * particleVelocityV(k);
                gridMass(gridx+i,gridy+j) = gridMass(gridx+i,gridy+j) + weight; % 这行似乎没用？
            end
        end
        
    end
    bound = 2;
    for i = 1:Nx
        for j = 1:Ny
            if gridMass(i,j) == 0
                gridVelocityU(i,j) = 0;
                gridVelocityV(i,j) = 0;
            else
                gridVelocityU(i,j) = gridVelocityU(i,j)/ gridMass(i,j);
                gridVelocityV(i,j) = gridVelocityV(i,j)/ gridMass(i,j);
            end
            gridVelocityV(i,j) = gridVelocityV(i,j) - dt * 0.98;
            if (i < bound) || (i > Nx - bound) || (j < bound) || (j > Ny - bound)
                gridVelocityU(i,j) = 0;
                gridVelocityV(i,j) = 0;
            end
        end
    end
    gridVelocityU(1,:) = 0;
    gridVelocityU(:,1) = 0;
    gridVelocityV(1,:) = 0;
    gridVelocityV(:,1) = 0;
    gridVelocityU(Nx,:) = 0;
    gridVelocityU(:,Ny) = 0;
    gridVelocityV(Nx,:) = 0;
    gridVelocityV(:,Ny) = 0;
    particleVelocityU = zeros(particleCount,1);
    particleVelocityV = zeros(particleCount,1);
    particleCx = zeros(particleCount,1);
    particleCy = zeros(particleCount,1);
    dt = 0.1;
    gridVelocityU(1,:) = 0.1;
    gridVelocityV(1,:) = -1;
    gridVelocityV(Nx,:) = 0;
    gridVelocityU(:,Ny) = 0;
    % Grid To Particle
    for k = 1:particleCount
        
        x = max(particlePosX(k) / dx + 1,1);
        y = max(particlePosY(k) / dx + 1,1);
        gridx = floor(x);
        gridy = floor(y);
        fx = x - gridx;
        fy = y - gridy;
        wx = [0.5 * (1.5 - fx).^2 0.75 - (fx - 1).^2 0.5*(fx - 0.5).^2];
        wy = [0.5 * (1.5 - fy).^2 0.75 - (fy - 1).^2 0.5*(fy - 0.5).^2];
        for j = 0:2
            for i = 0:2
                weight = wx(i+1) * wy(j+1);
                dposx = (i - fx) * dx;
                dposy = (j - fx) * dx;
                particleVelocityU(k) = particleVelocityU(k) + weight * gridVelocityU(gridx+i,gridy+j);
                particleVelocityV(k) = particleVelocityV(k) + weight * gridVelocityV(gridx+i,gridy+j);
                particleCx(k) = particleCx(k)  + 4 * weight * dposx * gridVelocityU(gridx+i,gridy+j) / dx / dx;
                particleCy(k) = particleCy(k)  + 4 * weight * dposy * gridVelocityV(gridx+i,gridy+j) / dx / dx;
            end
        end
        particlePosX(k) = particlePosX(k) + dt * particleVelocityU(k);
        particlePosY(k) = particlePosY(k) + dt * particleVelocityV(k);
    end
    
    time = time + 1;
    figure();
    xlim([0 1]);
    ylim([0,1]);
    pause(0.5);
    for i=1:length(particlePosX)
        xunit= particlePosX(i);
        yunit= particlePosY(i);
        hold on
        plot(xunit, yunit, 'Ob')
        hold off
    end
end

