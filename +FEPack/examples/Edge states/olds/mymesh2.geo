// Define the unit square
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(10) = {1, 2, 3, 4};
Plane Surface(11) = {10};

// Define a point at the center
Point(100) = {0.5, 0.5, 0, 0.1}; // The last value is a mesh size suggestion

// --- Distance field around center
Field[1] = Distance;
Field[1].NodesList = {100}; // Refine around this point

// --- Threshold field to control mesh size
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.02;  // Small elements near the point
Field[2].SizeMax = 0.2;   // Larger elements far from it
Field[2].DistMin = 0.1;   // Start refining within this radius
Field[2].DistMax = 0.3;   // Fully coarse beyond this distance

// Set it as the background mesh size field
Background Field = 2;
