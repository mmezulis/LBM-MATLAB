%% Nosacījumu uzdošana
global isObstacle Nx Ny dt dx X Y tau Uw XX YY;
Lx = 50; Nx = Lx; dx = Lx/Nx; X = linspace(dx/2, Lx-dx/2, Nx); %Režģa izmēri
Ly = 50; Ny = Ly/dx; Y = linspace(dx/2, Ly-dx/2, Ny);
[XX, YY] = meshgrid(X, Y);
endT = 1000; dt = 1; tt = 0:dt:endT; tau = 1;

% Šķērslis - cilinds apgabala centrā
isObstacle = (((XX' - Lx/2).^2 + (YY' - Ly/2).^2) < (Lx*0.2)^2) | (YY' < 0) | (YY' > Ly);

% Sākumnosacījumi
RHO = ones(Nx, Ny);
UX = zeros(Nx, Ny);
UY = zeros(Nx, Ny);

Uw = [0.001, 0]; % Ātrums uz augšējās un apakšējās sienas

f = initializeAlgorithm(RHO, UX, UY);
%% Cikls
T = struct('macro', 0, 'equilibrium', 0, 'collision', 0, 'streaming', 0);
for t=tt
    t
    tic
        % Calculate RHO and UX, UY
        [RHO, UX, UY, V] = calculate_macro(f);
    T.macro = T.macro + toc;

%         tic
%             visualize(RHO, UX, UY, V); 
%         T.visualize = toc;

    tic
        % calculate equilibrium
        fEq = equilibrium(UX, UY, RHO);
    T.equilibrium = T.equilibrium + toc;

    tic
        % collision step
        ftemp = collision(f, fEq);
    T.collision = T.collision + toc;

    tic
        %streaming (propagation) step
        f = streaming(ftemp, RHO);
    T.streaming = T.streaming + toc;
end

%% Straumēšanas solis
function f = streaming(ftemp, RHO)
    global b fIdx rhoMultiplier;

    f = ftemp(fIdx) + rhoMultiplier.*repmat(RHO, [1, 1, b]);
end

%% Makroskopisko lielumu aprēķināšana
function [RHO, UX, UY, U] = calculate_macro(f)
    global b c; 
    RHO = sum(f, 3);
    
    UX = zeros([size(f,1), size(f,2)]);
    UY = zeros([size(f,1), size(f,2)]);
    for ie=1:b
       UX(:, :) = UX(:, :) + f(:, :, ie)*c(ie, 1);
       UY(:, :) = UY(:, :) + f(:, :, ie)*c(ie, 2);
    end
    
    UX = UX./RHO;
    UY = UY./RHO;
    
    U = sqrt(UX.^2 + UY.^2);
end

%% Līdzsvara sadalījuma aprēķināšana
function fEq = equilibrium(UX, UY, RHO)
    global b c cs2inv cs4inv w Nx Ny;
    fEq = zeros(Nx, Ny, b);
    uu = UX.^2 + UY.^2;
    for ie = 1:b
        uci = UX*c(ie, 1) + UY*c(ie, 2);
        fEq(:, :, ie) = w(ie)*RHO.*(1 + uci*cs2inv + uci.^2*cs4inv/2 - uu*cs2inv/2);
    end
end

%% Sadursmju solis
function ftemp = collision(f, fEq)
    global Nx Ny b stepProp stepPropInv isObstacleFull;

    ftemp = zeros(Nx, Ny, b);
    for ie = 1:b
        ftemp(:, :, ie) = f(:, :, ie)*stepPropInv + fEq(:, :, ie)*stepProp;
    end
end

%% Logu sagatavošana izvaddatiem
function initOutput(RHO, UX, UY, U)
    global isObstacle I_RHO I_V I_CURL X Y;
    figure(1);
    I_RHO = imagesc(X, Y, RHO', 'AlphaData', ~isObstacle', [0.8, 1.2]);
    colormap(jet);
    colorbar;
    title('Density');
    figure(2);
    I_V = imagesc(X, Y, U', 'AlphaData', ~isObstacle', [0, 0.01]);
    colormap(jet);
    colorbar;
    title('Velocity');
    figure(3);
    I_CURL = imagesc(X, Y, curl(UX, UY)', 'AlphaData', ~isObstacle', [-0.0001, 0.0001]);
    colormap(jet);
    colorbar;
    title('Vorticity');
end

%% Izvaddatu logu atjaunināšana ar aktuālajiem datiem
function visualize(RHO, UX, UY, V)
    global I_RHO I_V I_CURL;
    set(I_RHO, 'CData', RHO');
    set(I_V, 'CData', V');
    set(I_CURL, 'CData', curl(UX, UY)');
    drawnow limitrate;
end

%% Pretējo virzienu 'iekešošana'
function precalculateInverseDirections()
    global cNorm b invDirections;
    invDirections = zeros(1, b);
    for ie=1:b
        invDirections(ie) = find(ismember(cNorm, -cNorm(ie, :), 'rows'), 1, 'first'); 
    end
end

%% Algoritma inicializācija pie dotajiem SN
function f = initializeAlgorithm(RHO, UX, UY)
    global cs dx dt stepProp stepPropInv tau cs2inv cs4inv cNorm c b w;
    cs = sqrt(dx^2/dt^2/3);
    w = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36;
    
    stepProp = dt/tau; stepPropInv = 1 - dt/tau;
    cs2inv = 1/cs^2;
    cs4inv = 1/cs^4;

    cNorm = [0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1];
    c = cNorm*dx/dt;
    b = length(c);
    
    precalculateInverseDirections();
    initializeStreamIndices();

    % Sākumā sadalījums ir līdzsvara stāvoklī
    f = equilibrium(UX, UY, RHO);
    initOutput(RHO, UX, UY, sqrt(UX.^2 + UY.^2));
end

function initializeStreamIndices()
    global isObstacle Nx Ny b c fIdx rhoMultiplier cNorm invDirections Uw cs w;
    fIdx = zeros(Nx, Ny, b);
    rhoMultiplier = zeros(Nx, Ny, b);
    % 'Ievelkam' sadalijumus no kaimiņu mezgliem
    s = size(fIdx);
    % Iterējam cauri katram režģa mezglam
    for ix=1:Nx
        for iy=1:Ny
            % Iterējam cauri katram virzienam, kas iziet no šūnas
            for ie=1:b
                if isObstacle(ix, iy)
                    % Ja ir šķērslis, tad atstājam esošo vērtību
                    fIdx(ix, iy, ie) = sub2ind(s, ix, iy, ie);
                    continue;
                end
                % Periodisks RN pa x
                from_ix = mod(ix - cNorm(ie, 1) - 1, Nx) + 1;
                from_iy = iy - cNorm(ie, 2);

                if (from_iy == 0) || (from_iy == Ny+1)
                    % Ja tiek ievilktas no sienas, tad ņemam no šīs pašas
                    % šūnas, tikai pretējā virzienā
                    fIdx(ix, iy, ie) = sub2ind(s, ix, iy, invDirections(ie));
                    
                    % Kustīgas sienas RN
                    rhoMultiplier(ix, iy, ie) = -2*w(ie)*dot(-c(ie,:), Uw)/cs^2;
                elseif isObstacle(from_ix, from_iy)
                    % Ja tiek ievilktas no šķēršļa, tad ņemam no šīs pašas
                    % šūnas, tikai pretējā virzienā
                    fIdx(ix, iy, ie) = sub2ind(s, ix, iy, invDirections(ie));
                else
                    % Ja ir no kurienes ievilkt
                    fIdx(ix, iy, ie) = sub2ind(s, from_ix, from_iy, ie);
                end
            end
        end
    end
end
