Point(1) = {0, 0, 0};
Point(2) = {1000, 0, 0};
Point(3) = {1000, 1000, 0};
Point(4) = {0, 1000, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Surface(7) = {6};
Physical Line(8) = {1};
Physical Line(9) = {3};
Physical Line(10) = {4};
Physical Line(11) = {2};

