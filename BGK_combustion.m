global isObstacle e invDirections Lx Ly b dt isObstacle_full omega cs cs2inv cs4inv X Y F_X F_Y boundaryNodeIndex;
Nt = 50000;
dx = 1; %Attalums starp shunam
dy = dx;
dt = 1;
cs = sqrt(dx^2/dt^2/3);
cs2inv = 1/cs^2;
cs4inv = 1/cs^4;

% Obstacles
% LOAD FROM FILE
isObstacle = 1-imread('porous_half_sm.bmp')';
Lx = size(isObstacle,1);
Ly = size(isObstacle,2);

isObstacle_full = padarray(isObstacle, [1, 1], 1);
X = linspace(-dx/2, Lx+dx/2, Lx+2);
Y = linspace(-dy/2, Ly+dy/2, Ly+2); %Rezgis ar sienam uz 0 un L
[XX, YY] = meshgrid(X, Y);

% Uzzimet normales vektorus:
% ===============
% [EdgeX, EdgeY] = edgeNormal(isObstacle_full);
% imagesc(X, Y, isObstacle_full);
% colormap(autumn);
% hold on;
% quiver(XX, YY, EdgeX, EdgeY);
% hold off;
% return;

prepareBoundaryNodeIndex(isObstacle_full);
D = 0.3; % SV koeficients
tauF = 2;
tauG = thermalRelaxationTime(D);
tauC1 = thermalRelaxationTime(0.1); % Abu vielu difuzijas koeficients
tauC2 = thermalRelaxationTime(0.1);

% Areja speka matrica
% F_X = 0*(YY < 60 & XX > 40 & YY < 60 & YY > 40)*0.005;
% F_Y = 0*(XX < 60 & XX > 40 & YY < 60 & YY > 40)*0.002;
% F_X = padarray(F_X, [1, 1], 0);
% F_Y = padarray(F_Y, [1, 1], 0);
F_X = zeros(Lx+2, Ly+2);
F_Y = zeros(Lx+2, Ly+2);

% viskozitate
v = dt*cs^2*(tauF - 1/2);

RHO = ones(Lx+2, Ly+2);
G = 0.2*ones(Lx+2, Ly+2);
C1 = zeros(Lx+2, Ly+2);
C2 = zeros(Lx+2, Ly+2);
VX = zeros(Lx+2, Ly+2);
VY = zeros(Lx+2, Ly+2);

RHO(isObstacle_full == 1) = NaN;
VX(isObstacle_full == 1) = 0;
VY(isObstacle_full == 1) = 0;
C1(1:(Lx/2), :) = 1;
C2((Lx/2):end, :) = 1;
% Pa vidu 'dzirkstele'
% G(Lx/2,Ly/2) = 1;

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
g = thermalEquilibrium(VX, VY, G);
c1 = thermalEquilibrium(VX, VY, C1);
c2 = thermalEquilibrium(VX, VY, C2);
initOutput(RHO);

T = {};
% katra laika solii

G_RES = G;
C1_RES = C1;
C2_RES = C2;

for i=1:Nt
    tFrame = tic;
        tic
            % Calculate RHO and VX, VY
            [RHO, G, C1, C2, VX, VY, V] = calculate_macro(f, g, c1, c2);
        T.macro = toc;
        
        if(mod(i, 1000) == 0)
            G_RES(:, :, end+1) = G;
            C1_RES(:, :, end+1) = C1;
            C2_RES(:, :, end+1) = C2;
        end

        tic
            visualize(RHO, G, C1, C2, VX, VY, V);
        T.visualize = toc;
        
        tic
            F = forcingTerm(VX, VY, tauF);
            
            % Degsanas process
            k = 0.1*exp(-0.5./G); % Areniusa likums
            % Pienemam, ka degsanas process nerada siltuma vai
            % koncentracijas inerci, tikai rada siltumu
            
            q = 0.6; % ipatnejais sadegsanas siltums
            reactions = k.*(C1.*C2.^2);
            F_G = thermalSource(q*reactions);
            F_C1 = thermalSource(-reactions);
            F_C2 = thermalSource(-2*reactions);
        T.forcing = toc;

        tic
            % calculate equilibrium
            fEq = equilibrium(VX, VY, RHO);
            gEq = thermalEquilibrium(VX, VY, G);
            c1Eq = thermalEquilibrium(VX, VY, C1);
            c2Eq = thermalEquilibrium(VX, VY, C2);
        T.equilibirum = toc;

        tic
            % collision step
            ftemp = collision(f, fEq, F, tauF);
            gtemp = collision(g, gEq, F_G, tauG);
            c1temp = collision(c1, c1Eq, F_C1, tauC1);
            c2temp = collision(c2, c2Eq, F_C2, tauC2);
        T.collision = toc;

        tic
            %streaming (propagation) step
            f = streaming(ftemp, fEq);
            g = thermalStreaming(gtemp, gEq, G);
            c1 = thermalStreaming(c1temp, c1Eq, C1);
            c2 = thermalStreaming(c2temp, c2Eq, C2);
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
        % nobidam atpakal
        toAdd = circshift(toAdd, -e(ie, :));
        f(:, :, invDirections(ie)) = f(:, :, invDirections(ie)) + toAdd;
    end
end

function g = thermalStreaming(gtemp, gEq, C)
    global invDirections isObstacle_full e b omega;
    
    gtemp(repmat(isObstacle_full == 1, [1,1,b])) = 0;
    g = stream(gtemp, e);
    
    % Atspogulojam atpakal no malam un skersliem
    for ie=1:b
        % vertibas, kuras iegaja ieksa "sienaas", ar Dirihle RN
        toAdd = isObstacle_full.*g(:, :, ie);
        % Neimana RN uz visam sienam
        % =============
        toAdd = -toAdd + 2*omega(ie)*C.*isObstacle_full;
        
        % nobidam atpakal
        toAdd = circshift(toAdd, -e(ie, :));
        
        g(:, :, invDirections(ie)) = g(:, :, invDirections(ie)) + toAdd;
    end
end

function [RHO, C, C1, C2, VX, VY, V] = calculate_macro(f, g, c1, c2)
    global b e isObstacle_full F_X F_Y dt; 
    RHO = sum(f, 3);
    C = sum(g, 3);
    C1 = sum(c1, 3);
    C2 = sum(c2, 3);
    
    VX = zeros([size(f,1), size(f,2)]);
    VY = zeros([size(f,1), size(f,2)]);
    for ie=1:b
       VX(:, :) = VX(:, :) + f(:, :, ie)*e(ie, 1);
       VY(:, :) = VY(:, :) + f(:, :, ie)*e(ie, 2);
    end
    VX = VX + F_X*dt/2;
    VY = VY + F_Y*dt/2;
    
    VX = VX./RHO;
    VY = VY./RHO;
    
    V = sqrt(VX.^2 + VY.^2);
    
    RHO(isObstacle_full == 1) = NaN;
    C(isObstacle_full == 1) = NaN;
    
    % Neimana RN
    C = setBoundaryConcentration(C);
    C1 = setBoundaryConcentration(C1);
    C2 = setBoundaryConcentration(C2);
    
    VX(isObstacle_full == 1) = 0;
    VY(isObstacle_full == 1) = 0;
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
        fEq(1, :, ie) = omega(ie)*RHO(1,:).*uci(1,:)*cs2inv;
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
    I1 = imshow(RHO', [0, 1]);
    set(I1, 'AlphaData', ~isObstacle_full');
    colormap(jet);
    colorbar;
    title('Temperature');
    figure(2);
    I2 = imshow(RHO', [0, 1]);
    set(I2, 'AlphaData', ~isObstacle_full');
    title('O2 Concentration');
    colormap(jet);
    colorbar;
    figure(3);
    I3 = imshow(RHO', [0, 1]);
    set(I3, 'AlphaData', ~isObstacle_full');
    title('H2 Concentration');
    colormap(jet);
    colorbar;
end

function visualize(RHO, T, C1, C2, VX, VY, V)
    global I1 I2 I3;
    set(I1, 'CData', T');
    set(I2, 'CData', C1');
    set(I3, 'CData', C2');
    drawnow limitrate;  
end

function F = forcingTerm(VX, VY, tau)
    global b e cs2inv cs4inv F_X F_Y omega;
    F = zeros([size(VX), b]);
    for ie = 1:b
        uci = VX*e(ie, 1) + VY*e(ie, 2);
        F(:, :, ie) = (1-1/tau/2)*omega(ie)*(((e(ie,1)- VX)*cs2inv + uci*e(ie,1)*cs4inv).*F_X + ((e(ie,2)- VY)*cs2inv + uci*e(ie,2)*cs4inv).*F_Y);
    end
end

function Q = thermalSource(q)
    global omega b;
    Q = zeros([size(q), b]);
    for ie=1:b
        Q(:,:,ie) = omega(ie)*q;
    end
end

function tauG = thermalRelaxationTime(D)
    global cs dt;
    syms tau;
    tauG = double(solve(D == cs^2*(tau - dt/2)));
end


% returns vectors pointing away from edges of obstacle matrix
function [NormY, NormX] = edgeNormal(I)
    [Gmag, Gdir] = imgradient(I, 'prewitt');
    DirX = sin(deg2rad(Gdir)); DirY = -cos(deg2rad(Gdir));
    DirX((DirX < 0) & (DirX > -1e-6)) = 0; DirY((DirY < 0) & (DirY > -1e-6)) = 0;
    NormX = sign(DirX.*Gmag).*I; NormY = sign(DirY.*Gmag).*I;
end

function prepareBoundaryNodeIndex(isObstacle_full)
    global boundaryNodeIndex;
    [NormY, NormX] = edgeNormal(isObstacle_full);
    NormX = NormX'; NormY = NormY';
    boundaryNodeIndex = zeros(size(isObstacle_full));
    for ix = 1:size(isObstacle_full, 1)
        for iy = 1:size(isObstacle_full, 2)
            boundaryNodeIndex(ix, iy) = sub2ind(size(isObstacle_full), ix+NormX(iy, ix), iy+NormY(iy, ix)); 
        end
    end
end

function C = setBoundaryConcentration(C0)
    global boundaryNodeIndex;
    C = C0(boundaryNodeIndex);
    % Dazos gadijumos normales vektors var noradit uz mezglu, kurs atrodas
    % ieksa skerslii. Tada gadijuma nepieciesams C v?rt?bu panemt no tas,
    % kas jau ir skersli ieksa. Atkariba no sadu kesizu garuma, so var
    % nakties atkartot vel dazas reizes.
    C = C(boundaryNodeIndex);
end