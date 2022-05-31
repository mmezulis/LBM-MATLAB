global isObstacle e invDirections Lx Ly b dt isObstacle_full omega cs cs2inv cs4inv X Y;
Nt = 50000;
dx = 1; %Attalums starp shunam
dy = dx;
dt = 1;
cs = sqrt(dx^2/dt^2/3);
cs2inv = 1/cs^2;
cs4inv = 1/cs^4;

% Obstacles

% DRAW BY HAND
% Lx = 400;
% Ly = 400;
% isObstacle = zeros(Lx, Ly);
% X = linspace(dx/2, Lx-dx/2, Lx);
% Y = linspace(dy/2, Ly-dy/2, Ly); %Rezgis ar sienam uz 0 un L
% 
% while true
%     h = figure(1);
%     imagesc(X, Y, isObstacle');
%     rect = getrect(h);
%     x1 = rect(1); y1 = rect(2); w=rect(3); h = rect(4);
%     x2 = x1 + w; y2 = y1 + h;
%     [d, x1] = min(abs(X-x1)); [d, x2] = min(abs(X-x2));
%     [d, y1] = min(abs(Y-y1)); [d, y2] = min(abs(Y-y2));
%     if(h==0 && w == 0)
%         break
%     end
%     isObstacle(x1:x2, y1:y2) = 1;
% end

% LOAD FROM FILE
isObstacle = gpuArray(imread('valves_xl.bmp')');
Lx = size(isObstacle,1);
Ly = size(isObstacle,2);

isObstacle_full = padarray(isObstacle, [1, 1], 1);

X = linspace(dx/2, Lx-dx/2, Lx);
Y = linspace(dy/2, Ly-dy/2, Ly); %Rezgis ar sienam uz 0 un L
[XX, YY] = meshgrid(X, Y);

tauF = 0.6;

% Areja speka matrica
F = gpuArray.zeros(Lx+2, Ly+2);

% viskozitate
v = dt*cs^2*(tauF - 1/2);

RHO = gpuArray.ones(Lx+2, Ly+2);
VX = gpuArray.zeros(Lx+2, Ly+2);
VY = gpuArray.zeros(Lx+2, Ly+2);

RHO(isObstacle_full == 1) = NaN;
VX(isObstacle_full == 1) = 0;
VY(isObstacle_full == 1) = 0;

e = gpuArray([0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1]);
b = length(e); %2DQ9

invDirectionIndex = @(ie) find(ismember(e, -e(ie, :), 'rows'), 1, 'first');
%precalculate inverse directions for each ie
invDirections = gpuArray.zeros(1, b);
for ie=1:b
    invDirections(ie) = invDirectionIndex(ie); 
end

omega = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36;

colormap jet;

% Sakuma atrodas lidzsvara stavokli
f = equilibrium(VX, VY, RHO);
initOutput(RHO);

T = {};
% katra laika solii
for i=2:Nt
    tFrame = tic;
        tic
            % Calculate RHO and VX, VY
            [RHO, VX, VY, V] = calculate_macro(f);
        T.macro = toc;

        tic
            visualize(RHO, VX, VY, V);
        T.visualize = toc;

        tic
            % calculate equilibrium
            fEq = equilibrium(VX, VY, RHO);
        T.equilibirum = toc;

        tic
            % collision step
            ftemp = collision(f, fEq, tauF);
        T.collision = toc;

        tic
            %streaming (propagation) step
            f = streaming(ftemp, fEq);
        T.streaming = toc;
    T.total = toc(tFrame)
end

function f_res = stream(f, e)
    global b;
    % streamo katru no virzieniem
    f_res = gpuArray.zeros(size(f));
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
        toAdd(1,2:end-1) = toAdd(1,2:end-1) - 2*fEq(1, 2:end-1, ie);
        toAdd(end,2:end-1) = toAdd(end,2:end-1) - 2*fEq(end, 2:end-1, ie);
        
        % nobidam atpakal
        toAdd = circshift(toAdd, -e(ie, :));
        
        f(:, :, invDirections(ie)) = f(:, :, invDirections(ie)) + toAdd;
    end
end

function [RHO, VX, VY, V] = calculate_macro(f)
    global b e isObstacle_full dt; 
    RHO = sum(f, 3);
    
    VX = gpuArray.zeros([size(f,1), size(f,2)]);
    VY = gpuArray.zeros([size(f,1), size(f,2)]);
    for ie=1:b
       VX(:, :) = VX(:, :) + f(:, :, ie)*e(ie, 1);
       VY(:, :) = VY(:, :) + f(:, :, ie)*e(ie, 2);
    end
    
    VX = VX./RHO;
    VY = VY./RHO;
    
    V = sqrt(VX.^2 + VY.^2);
    
    RHO(isObstacle_full == 1) = NaN;
    
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
    VX(1,2:end-1) = 0.01;
    VY(1,2:end-1) = 0;
    VX(end, 2:end-1) = 0.01;
    VY(end, 2:end-1) = 0;
    RHO(1, 2:end-1) = RHO(2, 2:end-1);
    RHO(end, 2:end-1) = RHO(end-1, 2:end-1);
end

function fEq = equilibrium(VX, VY, RHO)
    global b e cs2inv cs4inv omega;
    fEq = gpuArray.zeros([size(RHO), b]);
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
        fEq(end, :, ie) = omega(ie)*RHO(end,:).*uci(end,:)*cs2inv;
    end
end

function ftemp = collision(f, fEq, tau)
    global b isObstacle_full dt;
    
    %collision (relaxation) step
    ftemp = gpuArray.zeros(size(f));
    for ie = 1:b
        ftemp(2:end-1, 2:end-1, ie) = f(2:end-1, 2:end-1, ie)*(1-dt/tau) + fEq(2:end-1, 2:end-1, ie)*dt/tau;
    end
    ftemp(isObstacle_full == 1) = NaN;
end

function initOutput(RHO)
    global isObstacle_full I1 I2 I3;
    figure(1);
    I1 = imshow(RHO', [0, 0.3]);
    set(I1, 'AlphaData', gather(~isObstacle_full'));
    colormap(jet);
    colorbar;
    title('Velocity');
    figure(2);
    I2 = imshow(RHO', [-0.001, 0.001]);
    set(I2, 'AlphaData', gather(~isObstacle_full'));
    title('Vorticity');
    colormap(jet);
    colorbar;
    figure(3);
    I3 = imshow(RHO', [0, 2]);
    set(I3, 'AlphaData', gather(~isObstacle_full'));
    title('Density');
    colormap(jet);
    colorbar;
end

function visualize(RHO, VX, VY, V)
    global X Y isObstacle_full I1 I2 I3;
    set(I1, 'CData', gather(V'));
    set(I2, 'CData', gather(curl(VX, VY)'));
    set(I3, 'CData', gather(RHO'));
    drawnow limitrate;  
end