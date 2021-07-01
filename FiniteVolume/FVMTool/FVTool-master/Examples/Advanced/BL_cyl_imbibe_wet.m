% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
% Written by Ali A. Eftekhari
% Last checked: June 2021
% Spontaneous imbibition in 2D cylindrical coordinate
% the code works fine and model the effect of wettability
% it is slow for parameter estimation purposes, otherwise fine
% 
clc
%% define the geometry
Nx = 20; % number of cells in x direction
Ny = 50; % number of cells in y direction
W = 0.02; % [m] length of the domain in x direction
H = 0.15; % [m] length of the domain in y direction
% m = createMesh1D(Nx, W);
m = createMeshCylindrical2D(Nx, Ny, W, H); % creates a 2D mesh
%% define the physical parametrs
krw0_ww = 0.3;
krw0_ow = 1.0;
kro0_ww = 0.6;
kro0_ow = 0.76;
nw_ww = 2.4;
nw_ow= 2.4;
no_ww = 2.0;
no_ow= 2.0;
sor_ww=0.1;
sor_ow=0.12;
swc_ww=0.09;
swc_ow=0.09;
SF0=0.3;
SF=createFaceVariable(m, SF0); % 1 is water wet, 0 is oil wet
krw0=krw0_ww*SF+krw0_ow*(1-SF);
kro0=kro0_ww*SF+kro0_ow*(1-SF);
sor=sor_ww*SF+sor_ow*(1-SF);
swc=swc_ww*SF+swc_ow*(1-SF);
no= no_ww*SF+no_ow*(1-SF);
nw= nw_ww*SF+nw_ow*(1-SF);
sws=@(sw, sor, swc)((sw>swc).*(sw<1-sor).*(sw-swc)./(1-sor-swc)+(sw>=1-sor));
kro=@(sw, kro0, sor, swc, no)((sw>=swc).*kro0.*(1-sws(sw, sor, swc)).^no+(sw<swc).*(1+(kro0-1)./swc.*sw));
krw=@(sw, krw0, sor, swc, nw)((sw<=1-sor).*krw0.*sws(sw, sor, swc).^nw+(sw>1-sor).*(-(1-krw0)./sor.*(1.0-sw)+1.0));
dkrwdsw=@(sw, krw0, sor, swc, nw)((sw<=1-sor).*nw.*krw0.*(1./(1-sor-swc)).*sws(sw, sor, swc).^(nw-1)+(sw>1-sor).*((1-krw0)./sor));
dkrodsw=@(sw, kro0, sor, swc, no)((sw>=swc).*(-kro0.*no.*(1-sws(sw, sor, swc)).^(no-1))./(-swc-sor+1)+(sw<swc).*((kro0-1)./swc));
p0 = 100e5; % [bar] pressure
pin = 150e5; % [bar] injection pressure at the left boundary
u_in= 1.0/(24*3600); % [m/s] equal to 1 m/day
sw0=swc_ww+0.1;
% sw0(10:end-10, 10:end-10)=swc+0.2;
% sw0 = swc+0.1; % initial water saturation
sw_in = 1;
mu_oil = 2e-3; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 0.1e-12; % [m^2] average reservoir permeability
phi0 = 0.2; % average porosity

eps1=1e-6;
clx=0.05;
cly=0.05;
V_dp=0.01; % Dykstra-Parsons coef.
perm_val= k0;%field2d(Nx,Ny,k0,V_dp,clx,cly);
k=createCellVariable(m, perm_val);
phi=createCellVariable(m, phi0);
teta_ow=deg2rad(130);
teta_ww=deg2rad(20);
gama_ow=0.03; % N/m
labda=2.4; % for Ekofisk chalk
b=0.7;
r_ave=sqrt(k0/phi0); %  meter average pore diameter
pce0=2*gama_ow*cos(teta_ww)/r_ave; % Pa capillary entry pressure
% pc0=1.0e9;
teta=teta_ww*SF+teta_ow*(1-SF);
pce=createFaceVariable(m, pce0);
% sw0=swc+(1-labda*log(pc0/pce)+sqrt((-1+labda*log(pc0/pce))^2+4*swc/(1-swc)))/2*(1-swc);
dpc=@(sw, pce, swc, sor, teta)dpc_imb(sw, pce, swc, sor, teta, labda, b);
% sw_plot=linspace(0,1, 10000);
% plot(sw_plot, pc(sw_plot))
lw = geometricMean(k)/mu_water;
lo = geometricMean(k)/mu_oil;
% Lw = @(sw)(krw(sw));
% Lo = @(sw)(k/mu_oil*kro(sw));
% dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
% dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% left boundary pressure gradient
% BCp.left.a(:)=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)); BCp.left.b(:)=0; BCp.left.c(:)=-u_in;
% change the right boandary to constant pressure (Dirichlet)
% BCp.left.a(:)=0; BCp.left.b(:)=1; BCp.left.c(:)=p0;
BCp.right.a(:)=0; BCp.right.b(:)=1; BCp.right.c(:)=p0;
BCp.top.a(:)=0; BCp.top.b(:)=1; BCp.top.c(:)=p0;
BCp.bottom.a(:)=0; BCp.bottom.b(:)=1; BCp.bottom.c(:)=p0+100;
% change the left boundary to constant saturation (Dirichlet)
% BCs.left.a(:)=0; BCs.left.b(:)=1; BCs.left.c(:)=1.0-sor;
sw_pc0= sw_zero_pc_imb(swc_ww*SF0+swc_ow*(1-SF0), sor_ww*SF0+sor_ow*(1-SF0),...
    teta_ww*SF0+teta_ow*(1-SF0), labda, b);
BCs.right.a(:)=0; BCs.right.b(:)=1; BCs.right.c(:)=sw_pc0;
BCs.top.a(:)=0; BCs.top.b(:)=1; BCs.top.c(:)=sw_pc0;
BCs.bottom.a(:)=0; BCs.bottom.b(:)=1; BCs.bottom.c(:)=sw_pc0;
%% define the time step and solver properties
% dt = 1000; % [s] time step
% dt=(W/Nx)/u_in/20; % [s]
dt=1.0;
t_end = 1000*dt; % [s] final time
eps_p = 1e-7; % pressure accuracy
eps_sw = 1e-7; % saturation accuracy
%% define the variables
sw_old = createCellVariable(m, sw0, BCs);
p_old = createCellVariable(m, p0, BCp);
sw = sw_old;
oil_init=domainInt(1-sw_old);
p = p_old;
uw = -gradientTerm(p_old); % an estimation of the water velocity
%% start the main loop
% generate intial pressure profile (necessary to initialize the fully
% implicit solver)

t = 0;
dsw_alwd= 0.01;
rec_fact=0.0;

while (t<t_end)
% for i=1:5
    error_p = 1e5;
    error_sw = 1e5;
    % Implicit loop
%     while ((error_p>eps_p) || (error_sw>eps_sw))
    while(1)
        for i=1:3
            % calculate parameters
            pgrad = gradientTerm(p);
    %         pcgrad=gradientTerm(pc(sw));
            sw_face = upwindMean(sw, -pgrad); % average value of water saturation
            sw_grad=gradientTerm(sw);
            sw_ave=arithmeticMean(sw);
            pcgrad=funceval(dpc, sw_ave, pce, swc, sor, teta).*sw_grad;
            % solve for pressure at known Sw
            labdao = lo.*funceval(kro, sw_face, kro0, sor, swc, no);
            labdaw = lw.*funceval(krw, sw_face, krw0, sor, swc, nw);
    %         dlabdaodsw = lo.*funceval(dkrodsw, sw_face);
    %         dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
            labda = labdao+labdaw;
            % compute [Jacobian] matrices
            Mdiffp1 = diffusionTerm(-labda);
            RHSpc1=divergenceTerm(labdao.*pcgrad);
            [Mbcp, RHSbcp] = boundaryCondition(BCp);
            RHS1 = RHSpc1+RHSbcp; % with capillary
            p_new=solvePDE(m, Mdiffp1+Mbcp, RHS1);

            % solve for Sw
            pgrad = gradientTerm(p_new);
            uw=-labdaw.*pgrad;
            [Mbcsw, RHSbcsw] = boundaryCondition(BCs);
            RHS_sw=-divergenceTerm(uw);
            sw_new=solveExplicitPDE(sw_old, dt, RHS_sw, BCs, phi);

            error_p = max(abs((p_new.value(:)-p_old.value(:))./p_new.value(:)));
            error_sw = max(abs(sw_new.value(:)-sw.value(:)));
            
            p=p_new;
            sw=sw_new;
        end
        if error_sw>dsw_alwd
            dt=dt*(dsw_alwd/error_sw);
            p=p_old;
            sw=sw_old;
        else
            t=t+dt;
            p = p_new;
            sw = sw_new;
            p_old = p;
            sw_old = sw;
            dt=min(dt*(dsw_alwd/error_sw), 100*dt)
            break;
        end
    end
    rec_fact=[rec_fact, (oil_init-domainInt(1-sw))/oil_init];
    figure(1);visualizeCells(1-sw); drawnow;
end