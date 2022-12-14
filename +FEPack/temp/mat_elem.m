function [Mloc, Kloc, KthetaLoc, Gloc1, Gloc2] = mat_elem(S1, S2, S3, mu, rho)
% MAT_ELEM Matrices elementaires pour des EF de Lagrange P1
%
% SYNOPSIS [Mloc, Kloc, KthetaLoc] = mat_elem(S1, S2, S3, coeff, params)
%
% ENTREES * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                        (vecteurs reels 1x2)
%
%         * params : structure contenant les donnees associees au
%                    probleme
%
%         * coeffs : fonctions mu et rho associees au probleme
%                    (coeff.mu et coeff.rho sont des fonctions)
%
%         * opts   : Definit la maniere dont les matrices sont calculees
%                    'constant' si les coefficients sont constants et
%                    'variable' sinon
%
%
% SORTIES - Mloc  : matrice de masse elementaire (matrice 3x3)
%         - Kloc  : matrice de rigidite elementaire (matrice 3x3)
%         - KthetaLoc : matrice de rigidite directionnelle elementaire
%                       (matrice 3x3)
%         - Gloc1 : matrice associee a la forme \int (dx1(u) * v)
%         - Gloc2 : matrice associee a la forme \int (dx2(u) * v)
%
% NOTES (1) le calcul est effectue a l'aide d'une quadrature de
%          Gauss-Lobatto
% ========================================================================

if (nargin < 4)
  coeffs.mu  = @(x, y) 1.0;
else
  coeffs.mu = mu;
end

if (nargin < 5)
  coeffs.rho = @(x, y) 1.0;
else
  coeffs.rho = rho;
end


opts = 'variable'; % 'constant';


% Coordonnees des sommets du triangle
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);
if (abs(D) <= eps)
  error('l aire d un triangle est nulle!!!');
end

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = [y2 - y3, x3 - x2;...
        y3 - y1, x1 - x3;...
        y1 - y2, x2 - x1];

% Vecteur direction de coupe
theta = pi/3;
eT = [cos(theta); sin(theta)];

% Integrales
if strcmp(opts, 'constant')

    % Cas de coefficients constants
    Mloc = (abs(D)/24) * [2 1 1; 1 2 1; 1 1 2];
    Kloc = (0.5/abs(D)) * (norm * norm.');
    KthetaLoc = (0.5/abs(D)) * ((norm * eT) * (norm * eT).');

else

    % Cas de coefficients variables
    % Formule de quadrature de Gauss-Lobatto a 4 points
    P = [x2-x1 x3-x1; y2-y1 y3-y1] * [1/3 1/5 1/5 3/5; 1/3 1/5 3/5 1/5] + [x1; y1]*ones(1,4);

    % Matrices de rigidite elementaires
    AA_KK = (-9/32) * coeffs.mu(P(1,1), P(2,1)) + (25/96) * coeffs.mu(P(1,2), P(2,2)) + ...
            (25/96) * coeffs.mu(P(1,3), P(2,3)) + (25/96) * coeffs.mu(P(1,4), P(2,4));
    Kloc  = (1/abs(D)) * norm * AA_KK * norm.';
    KthetaLoc = (1/abs(D)) * (norm * eT) * AA_KK * (norm * eT).';

    % Matrice de masse elementaire
    w1 = @(x, y) (1/D) * ((y2-y3) * (x-x3) - (x2-x3) * (y-y3));
    w2 = @(x, y) (1/D) * ((y3-y1) * (x-x1) - (x3-x1) * (y-y1));
    w3 = @(x, y) (1/D) * ((y1-y2) * (x-x2) - (x1-x2) * (y-y2));

    A = @(x,y) [w1(x, y); w2(x, y); w3(x, y)] * coeffs.rho(x, y) * [w1(x, y), w2(x, y), w3(x, y)];
    Mloc = abs(D) * ((-9/32) * A(P(1,1), P(2,1)) + (25/96) * A(P(1,2), P(2,2)) + ...
                     (25/96) * A(P(1,3), P(2,3)) + (25/96) * A(P(1,4), P(2,4)));

end

Gloc1 = (1/6) * ones(3,1) * norm(:,1)';
Gloc2 = (1/6) * ones(3,1) * norm(:,2)';

end
