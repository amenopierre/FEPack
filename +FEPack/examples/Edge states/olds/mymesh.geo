h0 = 0.1;
is_structured = 0;
k1x =  3.627598728468436;
k1y =  6.283185307179586;
k2x =  3.627598728468436;
k2y = -6.283185307179586;

Kx = 0;
Ky = 4.188790204786391;

x1 = -Kx + k1x + 0.1; y1 = -Ky - 0.1; z1 = 0.0;
x2 = -Kx + k1x + 0.1; y2 =  Ky + 0.1; z2 = 0.0;
x3 =  Kx - k1x - 0.1; y3 =  Ky + 0.1; z3 = 0.0;
x4 =  Kx - k1x - 0.1; y4 = -Ky - 0.1; z4 = 0.0;

numNodesX = 16; numNodesY = 16;

domain_name = "rect";
side_name1 = "ymin";
side_name2 = "xmax";
side_name3 = "ymax";
side_name4 = "xmin";

// Points
P_1 = newp; Point(P_1) = {x1, y1, z1, h0};
P_2 = newp; Point(P_2) = {x2, y2, z2, h0};
P_3 = newp; Point(P_3) = {x3, y3, z3, h0};
P_4 = newp; Point(P_4) = {x4, y4, z4, h0};

domain_0 = {P_1[], P_2[], P_3[], P_4[]};

// Brillouin zone vertices
B_1 = newp; Point(B_1) = { Kx,        Ky,       0.0, h0};
B_2 = newp; Point(B_2) = {-Kx + k1x, -Ky + k1y, 0.0, h0};
B_3 = newp; Point(B_3) = { Kx + k2x,  Ky + k2y, 0.0, h0};
B_4 = newp; Point(B_4) = {-Kx,       -Ky,       0.0, h0};
B_5 = newp; Point(B_5) = { Kx - k1x,  Ky - k1y, 0.0, h0};
B_6 = newp; Point(B_6) = {-Kx - k2x, -Ky - k2y, 0.0, h0};



// Lines
L_1 = newl; Line(L_1) = {P_1, P_2}; Transfinite Line {L_1} = numNodesX;
L_2 = newl; Line(L_2) = {P_2, P_3}; Transfinite Line {L_2} = numNodesY;
L_3 = newl; Line(L_3) = {P_3, P_4}; Transfinite Line {L_3} = numNodesX;
L_4 = newl; Line(L_4) = {P_4, P_1}; Transfinite Line {L_4} = numNodesY;
domain_1 = {L_1[]};
domain_2 = {L_2[]};
domain_3 = {L_3[]};
domain_4 = {L_4[]};

BL_1 = newl; Line(BL_1) = {B_1, B_2}; Transfinite Line {BL_1} = numNodesX;
BL_2 = newl; Line(BL_2) = {B_2, B_3}; Transfinite Line {BL_2} = numNodesX;
BL_3 = newl; Line(BL_3) = {B_3, B_4}; Transfinite Line {BL_3} = numNodesX;
BL_4 = newl; Line(BL_4) = {B_4, B_5}; Transfinite Line {BL_4} = numNodesX;
BL_5 = newl; Line(BL_5) = {B_5, B_6}; Transfinite Line {BL_5} = numNodesX;
BL_6 = newl; Line(BL_6) = {B_6, B_1}; Transfinite Line {BL_6} = numNodesX;

// Loops
LL_1 = newll; Line Loop(LL_1) = {L_1, L_2, L_3, L_4};
LL_2 = newll; Line Loop(LL_2) = {BL_1, BL_2, BL_3, BL_4, BL_5, BL_6};


// Plane surfaces
S_1 = news; Plane Surface(S_1) = {LL_1};
// Line Loop(LL_2) In Surface(S_1);
Line {BL_1, BL_2, BL_3, BL_4, BL_5, BL_6} In Surface{S_1};

// Make the mesh structured if specified
If (is_structured == 1)
  Transfinite Surface {S_1};
EndIf

domain_5 = {S_1[]};

Physical Point(1) = domain_0[];
Physical Line("ymin") = domain_1[];
Physical Line("xmax") = domain_2[];
Physical Line("ymax") = domain_3[];
Physical Line("xmin") = domain_4[];
Physical Surface("rect")= domain_5[];

Mesh.Format = 50;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;