global isObstacle e invDirections Lx Ly b dt isObstacle_full omega cs cs2inv cs4inv X Y F_X F_Y;
Nt = 50000;
dx = 1; %Attalums starp shunam
dy = dx;
dt = 1;
cs = sqrt(dx^2/dt^2/3);
cs2inv = 1/cs^2;
cs4inv = 1/cs^4;

% Obstacles

% DRAW BY HAND
Lx = 100;
Ly = 100;
isObstacle = zeros(Lx, Ly); 
X = linspace(dx/2, Lx-dx/2, Lx);
Y = linspace(dy/2, Ly-dy/2, Ly); %Rezgis ar sienam uz 0 un L

isObstacle_full = padarray(isObstacle, [1, 1], 1);
[XX, YY] = meshgrid(X, Y);

D = 0.01; % SV koeficients
tauF = 1;
tauG = thermalRelaxationTime(D);


% Areja speka matrica
F_X = zeros(Lx+2, Ly+2);
F_Y = 0.0001 + zeros(Lx+2, Ly+2);

% Arejais siltuma efekts
G = zeros(Lx+2, Ly+2);

% viskozitate
v = dt*cs^2*(tauF - 1/2);

RHO = ones(Lx+2, Ly+2);
C = zeros(Lx+2, Ly+2);
VX = zeros(Lx+2, Ly+2);
VY = zeros(Lx+2, Ly+2);

RHO(isObstacle_full == 1) = NaN;
VX(isObstacle_full == 1) = 0;
VY(isObstacle_full == 1) = 0;
C((2*Lx/5):(3*Lx/5), (2*Ly/5):(3*Ly/5)) = 1;
% RHO((2*Lx/5):(3*Lx/5), (2*Ly/5):(3*Ly/5)) = 1.1;

% flowUpper = 0;
% flowLower = 0;

e = [0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1];
b = length(e); %2DQ9

invDirectionIndex = @(ie) find(ismember(e, -e(ie, :), 'rows'), 1, 'first');
%precalculate inverse directions for each ie
invDirections = zeros(1, b);
for ie=1:b
    invDirections(ie) = invDirectionIndex(ie); 
end

omega = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36;

colormap jet;

% Sakuma atrodas lidzsvara stavokli
f = equilibrium(VX, VY, RHO);
g = thermalEquilibrium(VX, VY, C)
initOutput(RHO);

T = {};

C_RES = [];
% katra laika solii
for i=1:Nt
    tFrame = tic;
        tic
            % Calculate RHO and VX, VY
            [RHO, C, VX, VY, V] = calculate_macro(f, g);
        T.macro = toc;
        
        if(mod(i,1000) == 1)
           C_RES(:, :, end+1) = C; 
        end

        tic
            visualize(RHO, C, VX, VY, V);
        T.visualize = toc;
        
        tic
            % Thermal expansion coefficient
            alpha = 3;
            % Reference density
            rho_0 = 1;
            % Reference temperature
            t_0 = 0;
            % Buoyancy force
            FB_X = -alpha*rho_0*(C-t_0).*F_X;
            FB_Y = -alpha*rho_0*(C-t_0).*F_Y;
            F = forcingTerm(F_X + FB_X, F_Y + FB_Y, VX, VY, tauF);
        T.forcing = toc;

        tic
            % calculate equilibrium
            fEq = equilibrium(VX, VY, RHO);
            gEq = thermalEquilibrium(VX, VY, C);
        T.equilibirum = toc;

        tic
            % collision step
            ftemp = collision(f, fEq, F, tauF);
            gtemp = collision(g, gEq, zeros(Lx+2, Ly+2, b), tauG);
        T.collision = toc;

        tic
            %streaming (propagation) step
            f = streaming(ftemp, fEq);
            g = thermalStreaming(gtemp, gEq);
        T.streaming = toc;
    T.total = toc(tFrame)
end

function f_res = stream(f, e)
    global b;
    % streamo katru no virzieniem
    f_res = zeros(size(f));
    for k = 1:b
        f_res(:, :, k) = circshift(f(:,:,k),e(k, :));
    end
end

function f = streaming(ftemp, fEq)
    global invDirections isObstacle_full e b;
    
    ftemp(isnan(ftemp)) = 0;
    f = stream(ftemp, e);
    
    % Atspogulojam atpakal no malam un skersliem
    for ie=1:b
        % vertibas, kuras iegaja ieksa "sienaas"
        toAdd = isObstacle_full.*f(:, :, ie);
        % Blivuma RN
%         toAdd(1,2:end-1) = -toAdd(1,2:end-1) + 2*fEq(1, 2:end-1, ie);
%         toAdd(end,2:end-1) = -toAdd(end,2:end-1) + 2*fEq(end, 2:end-1, ie);

        % Atruma RN
%         toAdd(1,2:end-1) = toAdd(1,2:end-1) - 2*fEq(1, 2:end-1, ie);
%         toAdd(end,2:end-1) = toAdd(end,2:end-1) - 2*fEq(end, 2:end-1, ie);
        
        % nobidam atpakal
        toAdd = circshift(toAdd, -e(ie, :));
        
        f(:, :, invDirections(ie)) = f(:, :, invDirections(ie)) + toAdd;
    end
end

function g = thermalStreaming(gtemp, gEq)
    global invDirections isObstacle_full e b;
    
    gtemp(repmat(isObstacle_full == 1, [1,1,b])) = 0;
    g = stream(gtemp, e);
    
    % Atspogulojam atpakal no malam un skersliem
    for ie=1:b
        % vertibas, kuras iegaja ieksa "sienaas", ar Dirihle RN
        toAdd = isObstacle_full.*g(:, :, ie);
        
        %Neimana RN
        toAdd(1, :) =  -g(1,:,ie) + 2*gEq(1, :, ie);
        toAdd(end, :) = -g(end,:,ie) + 2*gEq(end, :, ie);
        toAdd(:, 1) = -g(:,1,ie) + 2*gEq(:, 1, ie);
        toAdd(:, end) = -g(:,end,ie) + 2*gEq(:, end, ie);
        
        % Dirihle RN
%         toAdd(1,2:end-1) = -toAdd(1,2:end-1) + 2*gEq(1, 2:end-1, ie);
        % nobidam atpakal
        toAdd = circshift(toAdd, -e(ie, :));
        
        g(:, :, invDirections(ie)) = g(:, :, invDirections(ie)) + toAdd;
    end
end

function [RHO, C, VX, VY, V] = calculate_macro(f, g)
    global b e isObstacle_full F_X F_Y dt; 
    C = sum(g, 3);
    RHO = sum(f, 3);

    VX = zeros([size(f,1), size(f,2)]);
    VY = zeros([size(f,1), size(f,2)]);
    for ie=1:b
       VX(:, :) = VX(:, :) + f(:, :, ie)*e(ie, 1);
       VY(:, :) = VY(:, :) + f(:, :, ie)*e(ie, 2);
    end
    VX = VX + F_X*dt/2;
    VY = VY + F_Y*dt/2;
    
    VX = VX./(RHO);
    VY = VY./(RHO);
    
    V = sqrt(VX.^2 + VY.^2);
    
    RHO(isObstacle_full == 1) = NaN;
    C(isObstacle_full == 1) = 1;
    
    % Neimana RN
    C = padarray(C(2:end-1, 2:end-1),[1 1],'replicate');
    
    VX(isObstacle_full == 1) = 0;
    VY(isObstacle_full == 1) = 0;
    % Ar BLIVUMA RN
    % Ekstrapolejam atrumu
%     VX(1, 2:end-1) = 3*VX(2, 2:end-1)/2 - VX(3, 2:end-1)/2;
%     VY(1, 2:end-1) = 3*VY(2, 2:end-1)/2 - VY(3, 2:end-1)/2;
%     VX(end, 2:end-1) = 3*VX(end-1, 2:end-1)/2 - VX(end-2, 2:end-1)/2;
%     VY(end, 2:end-1) = 3*VY(end-1, 2:end-1)/2 - VY(end-2, 2:end-1)/2;
%     RHO(1, 2:end-1) = 1.5;
%     RHO(end, 2:end-1) = 1;

%     AR ATRUMA RN
%     VX(1,2:end-1) = 0;
%     VY(1,2:end-1) = -0.1;
%     VX(end, 2:end-1) = 0.01;
%     VY(end, 2:end-1) = 0;
%     RHO(1, 2:end-1) = RHO(2, 2:end-1);
    

%     RHO(end, 2:end-1) = RHO(end-1, 2:end-1);
end

function fEq = equilibrium(VX, VY, RHO)
    global b e cs2inv cs4inv omega;
    fEq = zeros([size(RHO), b]);
    uu = VX.^2 + VY.^2;
    for ie = 1:b
       uci = VX*e(ie, 1) + VY*e(ie, 2);
       fEq(:, :, ie) = omega(ie)*RHO.*(1 + uci*cs2inv + uci.^2*cs4inv/2 - uu*cs2inv/2);
       
       % Uz inlet/outlet ir savadak jarekina
       % BLIVUMA RN
%        fEq(1, :, ie) = omega(ie)*RHO(1,:).*(1 + uci(1,:).^2*cs4inv/2 - uu(1,:)*cs2inv/2);
%        fEq(end, :, ie) = omega(ie)*RHO(end,:).*(1 + uci(end,:).^2*cs4inv/2 - uu(end,:)*cs2inv/2);
        
%       ATRUMA RN
%         fEq(1, :, ie) = omega(ie)*RHO(1,:).*uci(1,:)*cs2inv;
%         fEq(end, :, ie) = omega(ie)*RHO(end,:).*uci(end,:)*cs2inv;
    end
end

function gEq = thermalEquilibrium(VX, VY, C)
    global b e cs2inv cs4inv omega;
    gEq = zeros([size(C), b]);
    uu = VX.^2 + VY.^2;
    for ie = 1:b
       uci = VX*e(ie, 1) + VY*e(ie, 2);
       gEq(:, :, ie) = omega(ie)*C.*(1 + uci*cs2inv + uci.^2*cs4inv/2 - uu*cs2inv/2);
       
       % Neimana uz visam sienam
       gEq(:, end, ie) = omega(ie)*C(:, end);
       gEq(:, 1, ie) = omega(ie)*C(:, end);
       gEq(1, :, ie) = omega(ie)*C(:, end);
       gEq(end, :, ie) = omega(ie)*C(:, end);

    end
end

function ftemp = collision(f, fEq, F, tau)
    global b isObstacle_full dt;
    
    %collision (relaxation) step
    ftemp = zeros(size(f));
    for ie = 1:b
        ftemp(2:end-1, 2:end-1, ie) = f(2:end-1, 2:end-1, ie)*(1-dt/tau) + fEq(2:end-1, 2:end-1, ie)*dt/tau;
    end
    %pieskaita arejo speku
    ftemp = ftemp + F*dt;
    ftemp(isObstacle_full == 1) = NaN;
end

function initOutput(RHO)
    global isObstacle_full I1 I2 I3;
    figure(1);
    I1 = imagesc(RHO', [0, 1]);
    set(I1, 'AlphaData', ~isObstacle_full');
    colormap(jet);
    colorbar;
    title('Concentration');
    figure(2);
    I2 = imagesc(RHO', [0, 0.02]);
    set(I2, 'AlphaData', ~isObstacle_full');
    title('Velocity');
    colormap(jet);
    colorbar;
    figure(3);
    I3 = imagesc(RHO', [0.8, 1.2]);
    set(I3, 'AlphaData', ~isObstacle_full');
    title('Density');
    colormap(jet);
    colorbar;
end

function visualize(RHO, C, VX, VY, V)
    global X Y isObstacle_full I1 I2 I3;
    set(I1, 'CData', C');
    set(I2, 'CData', V');   
    set(I3, 'CData', RHO');  
%     figure(4);
%     [XX, YY] = meshgrid([NaN, X, NaN], [NaN, Y, NaN]);
%     quiver(XX, YY, VX, VY);
    drawnow limitrate;  
end

function F = forcingTerm(F_X, F_Y, VX, VY, tau)
    global b e cs2inv cs4inv omega;
    F = zeros([size(VX), b]);
    for ie = 1:b
        uci = VX*e(ie, 1) + VY*e(ie, 2);
        F(:, :, ie) = (1-1/tau/2)*omega(ie)*(((e(ie,1)- VX)*cs2inv + uci*e(ie,1)*cs4inv).*F_X + ((e(ie,2)- VY)*cs2inv + uci*e(ie,2)*cs4inv).*F_Y);
    end
end

function tauG = thermalRelaxationTime(D)
    global cs dt;
    syms tau;
    tauG = double(solve(D == cs^2*(tau - dt/2)));
end

