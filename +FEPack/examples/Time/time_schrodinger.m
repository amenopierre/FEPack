% Solve 2D time dependent Schrödinger with Backward Euler
clear; clc;

HcObj = FEPack.applications.HoneycombObject('none', 'none');
HcObj.V = @(x) 10 * cos(x(:, 1:2) *  HcObj.dualVec1) +...
               10 * cos(x(:, 1:2) *  HcObj.dualVec2) +...
               10 * cos(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

HcObj.W = @(x)      sin(x(:, 1:2) *  HcObj.dualVec1) +...
                    sin(x(:, 1:2) *  HcObj.dualVec2) +...
                    sin(x(:, 1:2) * (HcObj.dualVec1  + HcObj.dualVec2));

HcObj.kappa = @(s) -1 + 2*FEPack.tools.domainwall(s+0.5);
delta = 1;

edge_param = 0.5*sqrt(2);
Qpot = @(x) zeros(size(x, 1), 1);
% Qpot = @(x) HcObj.V(x) + delta*HcObj.kappa(2*pi*delta * x(:, 1:2) * (-edge_param * HcObj.dualVec1 + HcObj.dualVec2)) .* HcObj.W(x);

initial_data = @(x) FEPack.tools.cutoff(sqrt(x(:, 1).^2 + x(:, 2).^2), -0.5, 0.5);

L = 20;
N = 256;
msh = FEPack.meshes.MeshRectangle(1, [-L L], [-L L], N, N);
dom = msh.domain('volumic');
edge0x = msh.domain('xmin');
edge1x = msh.domain('xmax');
edge0y = msh.domain('ymin');
edge1y = msh.domain('ymax');

% Finite Element matrices
u = FEPack.pdes.PDEObject; 
v = dual(u);

mat_gradu_gradv = FEPack.pdes.Form.intg(dom, grad2(u) * grad2(v));
mat_Qpot_u_v    = FEPack.pdes.Form.intg(dom, (Qpot * u) * v);
mat_u_v         = FEPack.pdes.Form.intg(dom, u * v);

ecs = ((((u|edge0x) - (u|edge1x)) == 0.0) &...
       (((u|edge0y) - (u|edge1y)) == 0.0));

ecs.applyEcs;
P = ecs.P;

% Time step
delta_t = 0.005;
t_final = 1;
t_current = 0;

%% Initialization
AA =  0.5i * delta_t * (mat_gradu_gradv + mat_Qpot_u_v) + mat_u_v;
BB = -0.5i * delta_t * (mat_gradu_gradv + mat_Qpot_u_v) + mat_u_v;
AA0 = P * AA * P';
BB0 = P * BB * P';
sol = initial_data(msh.points);

% %%
figure;
trisurf(msh.triangles, msh.points(:, 1), msh.points(:, 2), Qpot(msh.points));
shading interp;
view(2);

figure;
trisurf(msh.triangles, msh.points(:, 1), msh.points(:, 2), real(sol)); hold on;
% caxis([0, 1]);
shading interp;
view(2);
hLine = plot3([-L*(HcObj.vecPer1(1) + edge_param * HcObj.vecPer2(1)),...
                L*(HcObj.vecPer1(1) + edge_param * HcObj.vecPer2(1))],...
              [-L*(HcObj.vecPer1(2) + edge_param * HcObj.vecPer2(2)),...
                L*(HcObj.vecPer1(2) + edge_param * HcObj.vecPer2(2))],...
                [0 0], 'w-', 'LineWidth', 1.5);
axis([-L L -L L]);
set(gca, 'DataAspectRatio', [1 1 1]);
% ax = gca;
% ax.Children = [hLine; ax.Children(ax.Children ~= hLine)];
hold off;
pause;

while (true)

  rhs = BB0 * P * sol;
  sol = P' * (AA0 \ rhs);
  
  t_current = t_current + delta_t;

  trisurf(msh.triangles, msh.points(:, 1), msh.points(:, 2), real(sol)); hold on
  shading interp;
  hLine = plot3([-L*(HcObj.vecPer1(1) + edge_param * HcObj.vecPer2(1)),...
                  L*(HcObj.vecPer1(1) + edge_param * HcObj.vecPer2(1))],...
                [-L*(HcObj.vecPer1(2) + edge_param * HcObj.vecPer2(2)),...
                  L*(HcObj.vecPer1(2) + edge_param * HcObj.vecPer2(2))],...
                  [0 0], 'w-', 'LineWidth', 1.5);
  axis([-L L -L L]);
  set(gca, 'DataAspectRatio', [1 1 1]);
  view(2);
  title(['t = ', num2str(t_current, '%f')]);
  hold off;
  % caxis([0, 1]);

  if (t_current >= t_final)
    break;
  end

  pause(0.05);

end