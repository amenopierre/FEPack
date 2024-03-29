import FEPack.*

m = load('pregenMeshes/2D/unstruct_mesh_2D_10_positive.mat');

solbox = FEPack.applications.BVP.PerHalfGuideBVP(2, ...
  'semiInfiniteDirection', 1, ...
  'mesh', m.mesh, ...
  'volumeBilinearIntegrand', @(u, v) grad(u) * grad(v), ...
  'volumeLinearIntegrand', @(v) v, ...
  'boundaryConditions', @(u, xmin, ymin, ymax) ((u|xmin) == 0) && (((u|ymin) - (u|ymax)) == 0) && (((dn(u)|ymin) - (dn(u)|ymax)) == 0));