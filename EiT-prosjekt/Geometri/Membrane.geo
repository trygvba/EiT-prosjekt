cl1 = 0.4;
cl2 = 0.4;
Point(1) = {-1, -1, 10, cl1};
Point(2) = {-1, 1, 10, cl1};
Point(3) = {1, -1, 10, cl1};
Point(4) = {1, 1, 10, cl1};
Point(5) = {-1, -1, 4, cl2};
Point(6) = {-1, 1, 4, cl2};
Point(7) = {1, -1, 4, cl2};
Point(8) = {1, 1, 4, cl2};
Point(9) = {-1, -1, -4, cl2};
Point(10) = {-1, 1, -4, cl2};
Point(11) = {1, -1, -4, cl2};
Point(12) = {1, 1, -4, cl2};
Point(13) = {-1, -1, -10, cl1};
Point(14) = {-1, 1, -10, cl1};
Point(15) = {1, -1, -10, cl1};
Point(16) = {1, 1, -10, cl1};
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {2, 4};
Line(4) = {3, 4};
Line(5) = {5, 6};
Line(6) = {5, 7};
Line(7) = {6, 8};
Line(8) = {7, 8};
Line(9) = {9, 10};
Line(10) = {9, 11};
Line(11) = {10, 12};
Line(12) = {11, 12};
Line(13) = {13, 14};
Line(14) = {13, 15};
Line(15) = {14, 16};
Line(16) = {15, 16};
Line(17) = {1, 5};
Line(18) = {3, 7};
Line(19) = {4, 8};
Line(20) = {2, 6};
Line(21) = {8, 12};
Line(22) = {6, 10};
Line(23) = {5, 9};
Line(24) = {7, 11};
Line(25) = {10, 14};
Line(26) = {12, 16};
Line(27) = {11, 15};
Line(28) = {9, 13};
Line Loop(30) = {16, -15, -13, 14};
Plane Surface(30) = {30};
Line Loop(32) = {12, -11, -9, 10};
Plane Surface(32) = {32};
Line Loop(34) = {26, -15, -25, 11};
Plane Surface(34) = {34};
Line Loop(36) = {14, -27, -10, 28};
Plane Surface(36) = {36};
Line Loop(38) = {12, 26, -16, -27};
Plane Surface(38) = {38};
Line Loop(40) = {28, 13, -25, -9};
Plane Surface(40) = {40};
Line Loop(44) = {18, -6, -17, 2};
Plane Surface(44) = {44};
Line Loop(46) = {4, -3, -1, 2};
Plane Surface(46) = {46};
Line Loop(48) = {19, -8, -18, 4};
Plane Surface(48) = {48};
Line Loop(50) = {19, -7, -20, 3};
Plane Surface(50) = {50};
Line Loop(52) = {5, -20, -1, 17};
Plane Surface(52) = {52};
Line Loop(54) = {6, 8, -7, -5};
Plane Surface(54) = {54};
Line Loop(56) = {22, -9, -23, 5};
Plane Surface(56) = {56};
Line Loop(58) = {7, 21, -11, -22};
Plane Surface(58) = {58};
Line Loop(60) = {8, 21, -12, -24};
Plane Surface(60) = {60};
Line Loop(62) = {24, -10, -23, 6};
Plane Surface(62) = {62};
Surface Loop(63) = {44, 48, 50, 52, 46, 54}; // Top region
Volume(64) = {63};
Surface Loop(65) = {32, 54, 58, 60, 62, 56}; // Middle region
Volume(66) = {65};
Surface Loop(67) = {34, 38, 30, 40, 36, 32}; // Bottom region
Volume(68) = {67};
//Physical Groups:

Physical Surface(69) = {40, 36, 30, 38, 34, 58, 60, 62, 56, 50, 48, 44, 52, 46}; // Full boundary

Physical Volume(70) = {64}; // Region 1, top region
Physical Volume(71) = {66}; // Region 2, middle region
Physical Volume(72) = {68}; // Region 3, bottom region



