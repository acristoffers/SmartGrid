function F0 = findF0(A, C, F)
    % Finds F_0 with minimum-order such that Darouach's condition (4) in Ref. [1] is
    % satisfied for a triple (A,C,F_0). This is a MATLAB  implementation of
    % Algorithm 2 in Ref. [1]. See Ref. [1] for more details.
    %
    % Inputs:
    %   A        - dynamical matrix A of a dynamical system/network (size nxn)
    %   C        - output (measurement) matrix C (size qxn), which defines sensor
    %              nodes location
    %   F        - functional matrix F  (size rxn), which defines target nodes
    %              desired to be estimated
    %
    % Outputs:
    %   F0       - matrix F0 (size r0xn, where r0>=r), which defines the size
    %              minimum-order and structure of a functional observer
    %
    % References:
    %
    %   [1] A. N. Montanari, C. Duan, L. A. Aguirre, A. E. Motter. Functional
    %       observability and target state estimation in large-scale networks.

    % Copyright (C) 2020  Arthur Montanari
    %
    % This program is free software; you can redistribute it and/or modify it under
    % the terms of the GNU General Public License as published by the Free Software
    % Foundation; either version 2 of the License, or (at your option) any later
    % version.
    %
    % This program is distributed in the hope that it will be useful, but WITHOUT
    % ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    % FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    % details.
    %
    % The full text of the GNU General Public License can be found in the  file
    % LICENSE.

    % Last modified by Álan Crístoffer on 27/03/2021

    % Initialization
    F0 = zeros(size(A, 1), size(F, 2));
    F0(1:size(F, 1), :) = F;

    % Removes columns of A corresponding to sensor and target nodes.
    % See Note 1 for further explanation below.
    % A(:, [T; S]) = 0;

    % G is only composed by rows of CA, since C and F were removed (see Note 1
    % below).
    G = [C; C * A; F];
    M1 = [];

    % This is the while-loop in Algorithm 2 in Ref. [1]. This loop incrementally
    % augments F0 until the structural (generic) rank of [C; C*A; F0] is equal to
    % the structural rank of [C; C*A; F0; F0*A] -- satisfying Darouach's condition
    % (4).  Computes the structural rank of matrix G and [G; F0*A] using MATLAB's
    % function sprank. This is equivalent to building a corresponding bipartite
    % graph of a given matrix and finding the number of right-matched nodes via the
    % maximum matching algorithm.
    i = size(F, 1);
    while sprank(G) ~= sprank([G; F0(1:i, :) * A])
        Em = dmperm(G);
        M1 = unique([M1 find(Em)]);
        [~, M2] = find(F0(1:i, :) * A);
        Ca = setdiff(M2, M1);
        if isempty(Ca)
            [~, Ca] = find(F0);
            Ca = setdiff(find(Em), Ca);
        end
        x = Ca(1);
        i = i + 1;
        F0(i, x) = 1;
        G = [C; C * A; F0(1:i, :)];
    end

    F0 = F0(1:i, :);
end

%% Comments.

% Notation.
% Define G = [C; C*A; F0], and the bipartite graph B(V,X,E), where V is the set
% of nodes where each element corresponds to a row of G, X is a set of state
% nodes where each element also corresponds to a column of G, and (v_i,x_j) is
% an undirected edge in E if G(i,j)~=0.

% Note 1.
% If node j is a sensor or target node, then there is some row in G
% (corresponding to C or F0, respectively) such that the j-th column entry is a
% nonzero entry. Since, by assumption, this is a unique nonzero entry in this
% row, then the x_j \in X is a right-matched node, i.e., x_j \in M_1. Thus,
% sensor and target nodes are always right-matched nodes, i.e., it is not needed
% to check if they are right-matched or not.  Therefore, to reduce the
% dimensionality of the problem, we remove columns of A (and, thereby, of C*A,
% F0*A and G) associated with sensor and target nodes.
