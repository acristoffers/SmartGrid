%% Simulates attacks on a IEEE 118 powergrid.

% The MIT License (MIT);
%
% Copyright (c) 2021 Álan Crístoffer e Sousa
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% This script uses the IEEE 118 from [1]

disp("Initializing...");
model = pg_setup();

disp("Creating residual generators.");
model = pg_generators(model);

disp("Generating residuals on post-simulation attacks.");
model = pg_attack_post(model);

disp("Plotting graph.");
pg_plotgraph(model);
pg_plotsystem(model);

function [model] = pg_setup()
    version = 1;
    %
    addpath lib;
    addpath powergrid;
    %
    if isfile("powergrid/pg_saved.mat")
        load("powergrid/pg_saved", "model");
        if model.version == version
            return
        end
    end
    %
    model.version = version;
    %
    load("powergrid/IEEE118pg.mat", "*");
    load("powergrid/IEEE118pg.mat", "gamma"); % There is a gamma function, so we need to be explicit.
    %
    N            = 2 * ng + nl; % number of Kuramoto oscillators
    n            = 3 * ng + nl; % number of state variables
    linear.N     = N;
    linear.n     = n;
    linear.nl    = nl;
    linear.ng    = ng;
    linear.K     = K;
    model.linear = linear;
    % Simulates non-linear system to get x_{eq}
    disp("Simulating non-linear system.");
    t      = 0:1e-3:10;
    x0     = [theta0; zeros(ng, 1); theta0; zeros(nl, 1)];
    [~, x] = ode23tb(@(t, theta) pg_kuramoto(t, theta, ng, nl, wref, K, H, P, D, gamma, false), t, x0);
    x0     = x(end, :)';
    [~, x] = ode23tb(@(t, theta) pg_kuramoto(t, theta, ng, nl, wref, K, H, P, D, gamma, true), t, x0);
    % Places PMUs (sensors)
    disp("Placing PMUs (sensors).");
    [model] = pg_placePMUs(model);
    % Calculates linear system model
    disp("Calculating linear system.");
    A                               = zeros(n, n);
    A(1:ng, :)                      = [zeros(ng, ng), eye(ng), zeros(ng, ng + nl)]; % matrix block [11 12 13]
    A(ng + 1:2 * ng, ng + 1:2 * ng) = -D(1:ng) ./ (2 * H) .* eye(ng);               % matrix block [   22   ]
    thetaeq                         = [x0(1:ng, 1); x0(2 * ng + 1:end, 1)];
    aux_ij                          = wref * K .* cos(ones(N, 1) * thetaeq' - thetaeq * ones(1, N) + gamma) .* (ones(N, N) ./ ([2 * H; D(ng + 1:end)]));
    aux_ii                          = -wref * eye(N) .* sum(K .* cos(ones(N, 1) * thetaeq' - thetaeq * ones(1, N) + gamma)) .* (ones(N, N) ./ ([2 * H; D(ng + 1:end)]));
    A(ng + 1:end, 1:ng)             = aux_ij(:, 1:ng) + aux_ii(:, 1:ng);             % matrix block [21      ; 31      ];
    A(ng + 1:end, 2 * ng + 1:end)   = aux_ij(:, ng + 1:end) + aux_ii(:, ng + 1:end); % matrix block [      23;       33];
    B                               = zeros(size(A, 1), 1);
    C                               = full(sparse(1:length(model.linear.S), model.linear.S, 1, length(model.linear.S), size(A, 1)));
    D                               = zeros(size(C, 1), size(B, 2));
    model.linear.ss                 = ss(A, B, C, D);
    model.linear.x0                 = x0;
    model.linear.u0                 = zeros(size(B, 2), 1);
    simulation.t                    = t;
    simulation.x                    = x;
    simulation.u                    = zeros(length(t), size(B, 2));
    simulation.x0                   = x0;
    model.simulation                = simulation;
    %
    save("powergrid/pg_saved", "model", "-v7.3");
end

function [model] = pg_placePMUs(model)
    ng = model.linear.ng;
    nl = model.linear.nl;
    K  = model.linear.K;
    % PMUs are randomly placed on 30% of the generator terminals and load buses.
    N_PMU = ceil(0.3 * (ng + nl)); % number of PMUs
    PMU   = datasample(1:(ng + nl), N_PMU, "Replace", false);
    S     = PMU + 2 * ng; % set of sensor
    % Graph
    G = graph((K + K') / 2); % finds power grid (undirected) graph
    %
    model.linear.S = S;
    model.linear.G = G;
end

function pg_plotgraph(model)
    G  = model.linear.G;
    ng = model.linear.ng;
    N  = model.linear.N;
    S  = model.linear.S;
    %
    figure;
    hold on;
    %
    plot(0, nan, "*g");
    plot(0, nan, "*b");
    plot(0, nan, "*m");
    legend(["Generator nodes", "Generator terminals", "PMUs"]);
    %
    h            = plot(G, "Layout", "force", "Iterations", 40);
    h.EdgeColor  = [0.619, 0.619, 0.619];
    h.LineWidth  = 2;
    h.MarkerSize = 7;
    highlight(h, 1:1:ng, "NodeColor", "g");     % generator nodes
    highlight(h, ng + 1:1:N, "NodeColor", "b"); % generators terminals
    highlight(h, S - ng, "NodeColor", "m");     % PMUs
    set(gca, 'visible', 'off');
    fix2latex;
    %
    for i = 1:size(model.simulations_post, 1)
        if length(findobj("type", "figure")) > 100
            break
        end
        for j = 1:size(model.simulations_post, 2)
            if length(findobj("type", "figure")) > 100
                break
            end
            f = figure;
            plot(model.simulation.t, model.simulations_post{i, j}.r);
            title(sprintf("Attack %d on sensor %d", j, model.simulations_post{i, j}.attacked_state));
            xlabel("Time (s)");
            ylabel("Residuals");
            fix2latex;
            exportgraphics(f, sprintf("pg_imgs/s%da%d.eps", model.simulations_post{i, j}.attacked_state, j));
            for k = 1:length(model.simulations_post{i, j}.w)
                % if length(findobj("type", "figure")) > 100
                %     break
                % end
                % figure;
                % z = model.simulations_post{i, j}.z{k} + (model.generators{k}.F * model.linear.x0)';
                % plot(model.simulation.t, z);
                % title(strcat("Observer for attack ", num2str(j), " on sensor ", num2str(i)));
                % xlabel("Time (s)");
                % ylabel("States");
            end
        end
    end
end

function pg_plotsystem(model)
    t  = model.simulation.t;
    x  = model.simulation.x;
    ng = model.linear.ng;
    %
    figure;
    %
    subplot(121);
    plot(t, x(:, ng + 1:2 * ng));
    title("System frequencies");
    xlabel("Time (s)");
    ylabel("States");
    %
    subplot(122);
    plot(t, [x(:, 1:ng) x(:, 2 * ng + 1:end)]);
    title("System phases");
    xlabel("Time (s)");
    ylabel("States");
end

function [model] = pg_generators(model)
    if ~isfield(model, "generators")
        A = model.linear.ss.A;
        B = model.linear.ss.B;
        C = model.linear.ss.C;
        % One residual generator for each sensor.
        model.generators = cell(length(model.linear.S), 1);
        for i = 1:length(model.linear.S)
            disp("    Creating RG " + num2str(i) + "/" + num2str(length(model.linear.S)));
            % Attack vector (detectable one)
            L = full(sparse(model.linear.S(i), 1, 1, size(A, 1), 1));
            F = findF0(A, C, C);
            model.generators{i} = residual_generator(A, B, C, F, L, 1);
        end
        save("powergrid/pg_saved", "model", "-v7.3");
    end
end

function [model] = pg_attack_post(model)
    if isfield(model, "simulations_post")
        return
    end
    % Attacked nodes:
    L           = cellfun(@(g) find(g.L), model.generators);
    simulations = cell(length(L), 3);
    for i = 1:length(L)
        disp("    Attacks on sensor " + num2str(L(i)));
        l  = L(i);
        x0 = model.simulation.x0;
        tf = length(model.simulation.t);
        t0 = floor(2 / 5 * tf);
        tf = floor(4 / 5 * tf);
        % Attack 1: copy another state
        disp("        Attack 1");
        j                         = mod(i, length(L)) + 1;
        xa                        = model.simulation.x;
        xa(t0:tf, l)              = xa(t0:tf, L(j));
        xa                        = xa - x0';
        [r, w, z]                 = pg_calculate_residuals(model, xa);
        simulation.r              = sparse(r);
        simulation.z              = cellfun(@(c) sparse(c .* (c > 1e-4)), z, 'UniformOutput', false);
        simulation.w              = cellfun(@(c) sparse(c .* (c > 1e-4)), w, 'UniformOutput', false);
        simulation.attacked_state = l;
        simulations{i, 1}         = simulation;
        % Attack 2: constant offset
        disp("        Attack 2");
        xa                        = model.simulation.x;
        xa(t0:tf, l)              = xa(t0:tf, l) + 0.1;
        xa                        = xa - x0';
        [r, w, z]                 = pg_calculate_residuals(model, xa);
        simulation.r              = sparse(r);
        simulation.z              = cellfun(@(c) sparse(c .* (c > 1e-4)), z, 'UniformOutput', false);
        simulation.w              = cellfun(@(c) sparse(c .* (c > 1e-4)), w, 'UniformOutput', false);
        simulation.attacked_state = l;
        simulations{i, 2}         = simulation;
        % Attack 3: multiplicative
        disp("        Attack 3");
        xa                        = model.simulation.x;
        xa(t0:tf, l)              = xa(t0:tf, l) * 1.2;
        xa                        = xa - x0';
        [r, w, z]                 = pg_calculate_residuals(model, xa);
        simulation.r              = sparse(r);
        simulation.z              = cellfun(@(c) sparse(c .* (c > 1e-4)), z, 'UniformOutput', false);
        simulation.w              = cellfun(@(c) sparse(c .* (c > 1e-4)), w, 'UniformOutput', false);
        simulation.attacked_state = l;
        simulations{i, 3}         = simulation;
    end
    model.simulations_post = simulations;
    save("powergrid/pg_saved", "model", "-v7.3");
end

function [r, w, z] = pg_calculate_residuals(model, x)
    r          = cell(length(model.generators), 1);
    q          = cell(length(model.generators), 1);
    z          = cell(length(model.generators), 1);
    t          = model.simulation.t;
    u          = model.simulation.u' - model.linear.u0;
    generators = model.generators;
    for i = 1:length(r)
        g      = [generators{i}];
        y      = g.C * x';
        w      = zeros(size(g.N, 1), 1);
        [~, w] = ode23tb(@(t, w) pg_residual_observer(t, w, model, g, y, u), t, w);
        q{i}   = w;
        z{i}   = (g.D * w' + g.E * y)';
        r{i}   = (g.G * w' + g.M * y)';
    end
    r = cell2mat(r');
    w = q;
end

function [dw] = pg_residual_observer(t, w, model, generator, y, u)
    dt = model.simulation.t(2) - model.simulation.t(1);
    T  = length(model.simulation.t);
    k  = max(min(floor(t / dt), T), 2);
    y  = y(:, k - 1);
    u  = u(:, k - 1);
    %
    N  = generator.N;
    J  = generator.J;
    H  = generator.H;
    %
    dw = N * w + J * y + H * u;
end
