// define a len variable
len = 20.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {0.0, len, 0.0, 1.0};
Point(3) = {len, len, 0.0, 1.0};
Point(4) = {len, 0.0, 0.0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};

Line Loop(15) = {1,2,3,4};

Plane Surface(20) = {15};

Physical Surface(21) = {20};
