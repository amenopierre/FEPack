clear; clc;
N = 2*2.^(1:4);
alphas = 0.1.^(0:4);% 0.25*0.5.^(1:1);
import FEPack.*

u_exact = @(x) (x.^2 / 2) - (x.^3 / 3);
f = @(x, alpha) -1 + 2*x + alpha*u_exact(x);

for idJ = 1:length(alphas)
err = zeros(length(N), 1);

for idI = 1:length(N)

    mesh = meshes.MeshSegment('uniform', 0, 1, N(idI));

    u = FEPack.pdes.PDEObject;
    v = dual(u);
    
    MM = FEPack.pdes.Form.intg(mesh.domain('volumic'), u * v);
    KK = FEPack.pdes.Form.intg(mesh.domain('volumic'), grad(u, 1) * grad(v, 1));
    
    b = MM * f(mesh.points(:, 1), alphas(idJ));

    sol = (alphas(idJ) * MM + KK) \ b;
    solex = u_exact(mesh.points(:, 1));

    err(idI) = sqrt(((sol - solex)' * MM * (sol - solex)) / (solex' * MM * solex));
    
end

loglog(N, err, 'rx'); hold on;
p = polyfit(log(N), log(err), 1);
loglog(N, exp(p(2))*N.^(p(1)));
% pause;
end