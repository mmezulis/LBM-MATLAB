%% Nosacījumu uzdošana
global isObstacle Nx Ny dt dx X Y tau Uw XX YY;
Lx = 100; Nx = Lx; dx = Lx/Nx; X = linspace(dx/2, Lx-dx/2, Nx); %Režģa izmēri
Ly = 100; Ny = Ly + 2; Y = linspace(-dx/2, Ly+dx/2, Ny);
[XX, YY] = meshgrid(X, Y);
endT = 10000; dt = 1; tt = 0:dt:endT; tau = 1;

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
        [RHO, UX, UY, V] = calculate_macro(f);
    T.macro = T.macro + toc;

    visualize(RHO, UX, UY, V);

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
        f = streaming(ftemp, fEq);
    T.streaming = T.streaming + toc;
end

%% Straumēšanas solis
function f = streaming(ftemp, fEq)
    global invDirections isObstacle cNorm b isObstacleFull;
    
    function f_res = stream(f, cNorm)
        % streamo katru no virzieniem
        shifted = arrayfun(@(k) circshift(f(:,:,k),cNorm(k, :)),1:size(f,3),'uni',0);
        % sakombine atpakal par 3d matricu
        f_res = cat(3, shifted{:});
    end
    ftemp(isObstacleFull) = 0;
    f = stream(ftemp, cNorm);
    
    % Atspogulojam atpakal no malam un skersliem
    for ie=1:b
        % vertibas, kuras iegaja ieksa "sienaas"
        toAdd = isObstacle.*f(:, :, ie);

        % Atruma RN
        % Kustiga apakseja un augseja siena
        toAdd(:, 1) = toAdd(:, 1) - 2*fEq(:, 1, ie);
        toAdd(:, end) = toAdd(:, end) - 2*fEq(:, end, ie);
        
        % nobidam atpakal
        toAdd = circshift(toAdd, -cNorm(ie, :));
        
        f(:, :, invDirections(ie)) = f(:, :, invDirections(ie)) + toAdd;
    end
end

%% Makroskopisko lielumu aprēķināšana
function [RHO, UX, UY, V] = calculate_macro(f)
    global b c isObstacle Uw; 
    RHO = sum(f, 3);
    
    UX = zeros([size(f,1), size(f,2)]);
    UY = zeros([size(f,1), size(f,2)]);
    for ie=1:b
       UX(:, :) = UX(:, :) + f(:, :, ie)*c(ie, 1);
       UY(:, :) = UY(:, :) + f(:, :, ie)*c(ie, 2);
    end
    
    UX = UX./RHO;
    UY = UY./RHO;

%     AR ATRUMA RN
    UX(:,1) = Uw(1);
    UY(:,1) = Uw(2);
    RHO(:, 1) = RHO(:, 2);
    UX(:, end) = Uw(1);
    UY(:, end) = Uw(2);
    RHO(:, end) = RHO(:, end-1);
    
    V = sqrt(UX.^2 + UY.^2);
end

%% Līdzsvara sadalījuma aprēķināšana
function fEq = equilibrium(UX, UY, RHO)
    global b c cs2inv cs4inv w Nx Ny;
    fEq = zeros(Nx, Ny, b);
    uu = UX.^2 + UY.^2;
    for ie = 1:b
        uci = UX*c(ie, 1) + UY*c(ie, 2);
        fEq(:, :, ie) = w(ie)*RHO.*(1 + uci*cs2inv + uci.^2*cs4inv/2 - uu*cs2inv/2);
       
        % ĀTRUMA RN
        fEq(:, 1, ie) = w(ie)*RHO(:,1).*uci(:,1)*cs2inv;
        fEq(:, end, ie) = w(ie)*RHO(:,end).*uci(:,end)*cs2inv;
    end
end

%% Sadursmju solis
function ftemp = collision(f, fEq)
    global Nx Ny b stepProp stepPropInv isObstacleFull;

    ftemp = zeros(Nx, Ny, b);
    for ie = 1:b
        ftemp(:, :, ie) = f(:, :, ie)*stepPropInv + fEq(:, :, ie)*stepProp;
    end
    ftemp(isObstacleFull) = NaN;
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
    global cs dx dt stepProp stepPropInv tau cs2inv cs4inv cNorm c b w isObstacleFull isObstacle;
    cs = sqrt(dx^2/dt^2/3);
    
    stepProp = dt/tau; stepPropInv = 1 - dt/tau;
    cs2inv = 1/cs^2;
    cs4inv = 1/cs^4;

    cNorm = [0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1];
    c = cNorm*dx/dt;
    b = length(c);
    
    isObstacleFull = repmat(isObstacle == 1, [1, 1, b]);
    precalculateInverseDirections();

    w = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36;

    % Sākumā sadalījums ir līdzsvara stāvoklī
    f = equilibrium(UX, UY, RHO);
    initOutput(RHO, UX, UY, sqrt(UX.^2 + UY.^2));
end
