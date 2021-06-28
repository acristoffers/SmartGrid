function [generator] = residual_generator(A, B, C, F, Li, l)
    % Creates a residual generator for the functional system (A, B, C, F) and
    % fault Li. Li is the matrix that maps the fault to x, probably a column
    % vector.

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

    cvx_begin sdp quiet
        dim = size(F, 1);
        %
        variable P(dim, dim) symmetric semidefinite
        variable Y(dim, size(C, 1))
        variable K(dim, size(C, 1))
        %
        U = F * Li * pinv(C * Li);
        V = eye(size(Y, 2), size(U, 2)) - (C * Li) * pinv(C * Li);
        %
        Ah = A * pinv(F);
        Ch = C * pinv(F);
        Eh = P * U + Y * V;
        %
        X = Ah' * F' * P - Ah' * C' * Eh' - Ch' * K' + P * F * Ah - Eh * C * Ah - K * Ch - l * eye(dim, size(Ah, 2));
        X = triu(X, 1) + triu(X, 1)' + diag(diag(X));
        W = -sqrt(l) * (P * F - Eh * C);
        %
        minimize(1)
        subject to
            [X W; W' -eye(size(W, 2))] <= 0;
            P >= 1e-6 * eye(size(P, 1));
    cvx_end
    K = value(P) \ value(K);
    Y = value(P) \ value(Y);
    E = U + Y * V;
    R = F - E * C;
    N = (R * A - K * C) * pinv(F);
    J = K + N * E;
    H = R * B;
    D = eye(dim);

    cvx_begin quiet
        variable Q(1, size(C,1))
        variable l
        M = Q * (eye(size(C, 1)) - C * pinv(F) * E);
        G = -Q * C * pinv(F) * D;
        maximize(l)
        subject to
            G * R + M * C == 0;
            M >= l * ones(size(M, 1), size(M, 2));
    cvx_end
    Q = Q ./ norm(M);
    M = Q * (eye(size(C, 1)) - C * pinv(F) * E);
    G = -Q * C * pinv(F) * D;

    generator.F = F;
    generator.C = C;
    generator.L = Li;
    generator.N = N;
    generator.J = J;
    generator.H = H;
    generator.D = D;
    generator.E = E;
    generator.G = G;
    generator.M = M;
end
