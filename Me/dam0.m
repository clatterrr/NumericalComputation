Nx = 32;
Ny = 32;
dx = 1 / Nx;

gridDensity = zeros(Nx,Ny);
gridVelocityU = zeros(Nx,Ny);
gridVelocityV = zeros(Nx,Ny);
gridDivergence = zeros(Nx,Ny);
gridDivergence2 = zeros(Nx,Ny);
gridPressure = zeros(Nx,Ny);

particleCount = 60;
particleVelocityU = zeros(particleCount,1);
particleVelocityV = zeros(particleCount,1);
particlePosX = zeros(particleCount,1);
particlePosY = zeros(particleCount,1);

for i = 1:particleCount
    particlePosX(i) = floor((i-1) / 10) * dx;
    particlePosY(i) = mod(i-1,10) * dx;
    particleVelocityU(i) = (0.5 - particlePosX(i)) / 10;
    particleVelocityV(i) = -1;
end
for i = 1:Nx/2
    for j = 1:Ny/2
        gridPressure(i,j) = 0.5 - j * dx;
    end
end
time = 0;
timeFinal = 100;
while(time < timeFinal)
    gridVelocityU = zeros(Nx,Ny);
    gridVelocityV = zeros(Nx,Ny);
    % particle To Grid Velocity
    for k = 1:particleCount
        
        x = max(particlePosX(k) / dx + 1,1);
        y = max(particlePosY(k) / dx + 1,1);
        gridx = floor(x);
        gridy = floor(y);
        gridVelocityU(gridx,gridy) = gridVelocityU(gridx,gridy) + particleVelocityU(k) / 4;
        gridVelocityV(gridx,gridy) = gridVelocityV(gridx,gridy) + particleVelocityV(k) / 4;
        
        gridVelocityU(gridx,gridy+1) = gridVelocityU(gridx,gridy+1) + particleVelocityU(k) / 4;
        gridVelocityV(gridx,gridy+1) = gridVelocityV(gridx,gridy+1) + particleVelocityV(k) / 4;
        
        gridVelocityU(gridx+1,gridy) = gridVelocityU(gridx+1,gridy) + particleVelocityU(k) / 4;
        gridVelocityV(gridx+1,gridy) = gridVelocityV(gridx+1,gridy) + particleVelocityV(k) / 4;
        
        gridVelocityU(gridx+1,gridy+1) = gridVelocityU(gridx+1,gridy+1) + particleVelocityU(k) / 4;
        gridVelocityV(gridx+1,gridy+1) = gridVelocityV(gridx+1,gridy+1) + particleVelocityV(k) / 4;
    end
    gridVelocityU(1,:) = 0;
    gridVelocityU(:,1) = 0;
    gridVelocityV(1,:) = 0;
    gridVelocityV(:,1) = 0;
    gridVelocityU(Nx,:) = 0;
    gridVelocityU(:,Ny) = 0;
    gridVelocityV(Nx,:) = 0;
    gridVelocityV(:,Ny) = 0;
    totalD = 0;
    for i = 2:Nx-1
        for j = 2:Ny-1
            gridDivergence(i,j) = ( gridVelocityU(i,j) -  gridVelocityU(i-1,j) + ...
                gridVelocityV(i,j) -  gridVelocityV(i,j-1))/dx ;
            totalD = totalD + gridDivergence(i,j);
        end
    end
    for k = 1:10000
        for i = 2:Nx-1
            for j = 2:Ny-1
                div = 2 * (dx*dx + dx*dx);
                gridPressure(i,j) = ((gridPressure(i-1,j) + gridPressure(i+1,j))*dx*dx + (gridPressure(i,j-1) + ...
                    gridPressure(i,j+1))*dx*dx - gridDivergence(i,j)* dx * dx * dx * dx)/ div;
            end
        end
        gridPressure(1,:) = gridPressure(2,:);
        gridPressure(:,1) = gridPressure(:,2);
        gridPressure(Nx,:) = gridPressure(Nx-1,:);
        gridPressure(:,Ny) = gridPressure(:,Ny-1);
    end
    
    for i = 1:Nx-1
        for j = 1:Ny-1
            gridVelocityU(i,j) = gridVelocityU(i,j) + (gridPressure(i,j) - gridPressure(i+1,j))/dx;
            gridVelocityV(i,j) = gridVelocityV(i,j) + (gridPressure(i,j) - gridPressure(i,j+1))/dx;
        end
    end
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            gridDivergence2(i,j) = ( gridVelocityU(i,j) -  gridVelocityU(i-1,j) + gridVelocityV(i,j) -  gridVelocityV(i,j-1))/dx ;
        end
    end
    particleVelocityU = zeros(particleCount,1);
    particleVelocityV = zeros(particleCount,1);
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
        v00 = gridVelocityU(gridx,gridy);
        v10 = gridVelocityU(gridx+1,gridy);
        v01 = gridVelocityU(gridx,gridy+1);
        v11 = gridVelocityU(gridx+1,gridy+1);
        particleVelocityU(k) =  gridVelocityU(gridx,gridy) * (1 - fx) * (1 - fy) + ...
            gridVelocityU(gridx+1,gridy) * fx * (1 - fy) + ...
            gridVelocityU(gridx,gridy+1) * (1 - fx) * fy + ...
            gridVelocityU(gridx+1,gridy+1) * fx * fy;
        
        v00 = gridVelocityV(gridx,gridy);
        v10 = gridVelocityV(gridx+1,gridy);
        v01 = gridVelocityV(gridx,gridy+1);
        v11 = gridVelocityV(gridx+1,gridy+1);
        particleVelocityV(k) =  gridVelocityV(gridx,gridy) * (1 - fx) * (1 - fy) + ...
            gridVelocityV(gridx+1,gridy) * fx * (1 - fy) + ...
            gridVelocityV(gridx,gridy+1) * (1 - fx) * fy + ...
            gridVelocityV(gridx+1,gridy+1) * fx * fy;
        
        particlePosX(k) = particlePosX(k) + dt * particleVelocityU(k);
        particlePosY(k) = particlePosY(k) + dt * particleVelocityV(k);
        
    end
    
    time = time + 1;
    if mod(time,10) ~= 0
        continue
    end
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

