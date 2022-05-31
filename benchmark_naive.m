global c cNorm isObstacle invDirections Nx Ny b cs dt dx w X Y tau Uw XX YY X Y;
Lx = 50; Nx = Lx; dx = Lx/Nx; X = linspace(dx/2, Lx-dx/2, Nx); %Režģa izmēri
Ly = 50; Ny = Ly; Y = linspace(dx/2, Ly-dx/2, Ny);
[XX, YY] = meshgrid(X, Y);
endT = 1000; dt = 1; tt = 0:dt:endT;
isObstacle = (XX' - Lx/2).^2 + (YY' - Ly/2).^2 < (Lx*0.2)^2;

cs = sqrt(dx^2/dt^2/3);
tau = 1;

RHO = ones(Nx, Ny);
UX = zeros(Nx, Ny);
UY = zeros(Nx, Ny);

cNorm = [0, 0; 1, 0; 0, 1; -1, 0; 0, -1; 1, 1; -1, 1; -1, -1; 1, -1]; %2DQ9
c = cNorm*dx/dt;
b = length(c);
w = [16, 4, 4, 4, 4, 1, 1, 1, 1]/36; % Svaru koeficienti

Uw = [0.001, 0]; % Atrums uz augsejas un apaksejas sienas

% Aprēķinam pretējos ātruma virzienus
invDirections = @(ie) find(ismember(cNorm, -cNorm(ie, :), 'rows'), 1, 'first');

% Sākumā atrodas līdzsvara stāvoklī
f = equilibrium(UX, UY, RHO);
initOutput(RHO, UX, UY, sqrt(UX.^2 + UY.^2));

T = struct('macro', 0, 'equilibrium', 0, 'collision', 0, 'streaming', 0);
for t=tt
    t
    tic
        % Calculate RHO and UX, UY
        [RHO, UX, UY, U] = calculate_macro(f);
    T.macro = T.macro + toc;

%     visualize(RHO, UX, UY, U); 

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
streamslice(X, Y, UX', UY'); % Uzzīmējam plūsmas līnijas

% Straumēšana
function f = streaming(ftemp, RHO)
    global isObstacle invDirections c cNorm b Nx Ny Uw w cs;
    f = zeros(size(ftemp));
    % Iterējam cauri katram režģa mezglam
    for ix=1:Nx
        for iy=1:Ny
            if isObstacle(ix, iy)
                continue;
            end
            % Iterējam cauri katram virzienam, kas iziet no šūnas
            for ie=1:b
                % Periodisks RN pa x
                new_ix = mod(ix + cNorm(ie, 1) - 1, Nx) + 1;
                new_iy = iy+cNorm(ie, 2);
                
                if (new_iy == 0) || (new_iy == Ny+1)
                    % Ja notiek sadursme ar sienu, tad atstarojam daļiņas
                    % atpakaļ, ņemot vērā sienas ātrumu
                    f(ix, iy, invDirections(ie)) = ftemp(ix, iy, ie) - 2*w(ie)*RHO(ix, iy)*dot(c(ie, :), Uw)/cs^2;
                elseif isObstacle(new_ix, new_iy)
                    % Ja notiek sadursme ar šķērsli, tad atstarojam daļiņas atpakaļ
                    f(ix, iy, invDirections(ie)) = ftemp(ix, iy, ie);
                else
                    % Ja nav sadursme ar sienu, tad straumējam
                    f(new_ix, new_iy, ie) = ftemp(ix, iy, ie);
                end
            end
        end
    end
end

function [RHO, UX, UY, U] = calculate_macro(f)
    global c Nx Ny; 
    RHO = zeros(Nx, Ny);
    UX = zeros(Nx, Ny);
    UY = zeros(Nx, Ny);
    U = zeros(Nx, Ny);
    
    mdot = @(a, b) dot(a(:), b(:)); % Palīgfunkcija darbam ar n-dimensionālu vektoru skalāro reizinājumu
    for ix = 1:Nx
        for iy = 1:Ny
            RHO(ix, iy) = sum(f(ix, iy, :));
            UX(ix, iy) = mdot(c(:, 1), f(ix, iy, :))/RHO(ix, iy);
            UY(ix, iy) = mdot(c(:, 2), f(ix, iy, :))/RHO(ix, iy);
            U(ix, iy) = sqrt(UX(ix, iy)^2 + UY(ix, iy)^2);
        end
    end
end

function fEq = equilibrium(UX, UY, RHO)
    global b c cs w Nx Ny isObstacle;
    fEq = zeros(Nx, Ny, b);
    for ix = 1:Nx
        for iy = 1:Ny
            if isObstacle(ix, iy)
               continue; 
            end
            u = [UX(ix, iy), UY(ix, iy)];
            for ie = 1:b
                ci = c(ie, :);
                fEq(ix, iy, ie) = w(ie)*RHO(ix, iy)*(1 + dot(u, ci)/cs^2 + dot(u, ci)^2/cs^4/2 - dot(u,u)/cs^2/2);
            end
        end
    end
end

% SADURSMJU SOLIS
function ftemp = collision(f, fEq)
    global isObstacle Nx Ny b tau dt;

    ftemp = zeros(Nx, Ny, b);
    for ix = 1:Nx
        for iy = 1:Ny
            if isObstacle(ix, iy)
               continue; 
            end
            for ie = 1:b
                ftemp(ix, iy, ie) = f(ix, iy, ie)*(1-dt/tau) + fEq(ix, iy, ie)*dt/tau;
            end
        end
    end
end

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

function visualize(RHO, UX, UY, U)
    global I_RHO I_V I_CURL;
    set(I_RHO, 'CData', RHO');
    set(I_V, 'CData', U');
    set(I_CURL, 'CData', curl(UX, UY)');
    drawnow limitrate;
end
