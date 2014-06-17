cl1 = 0.2;
membraneratio = 0.5;
Point(1) = {0, 0, 1, cl1};
Point(2) = {0, 1, 1, cl1};
Point(3) = {1, 0, 1, cl1};
Point(4) = {1, 1, 1, cl1};

Point(5) = {0, 0, 0.1, membraneratio*cl1};
Point(6) = {0, 1, 0.1, membraneratio*cl1};
Point(7) = {1, 0, 0.1, membraneratio*cl1};
Point(8) = {1, 1, 0.1, membraneratio*cl1};

Point(9) = {0, 0, -0.1, membraneratio*cl1};
Point(10) = {0, 1, -0.1, membraneratio*cl1};
Point(11) = {1, 0, -0.1, membraneratio*cl1};
Point(12) = {1, 1, -0.1, membraneratio*cl1};

Point(13) = {0, 0, -1, cl1};
Point(14) = {0, 1, -1, cl1};
Point(15) = {1, 0, -1, cl1};
Point(16) = {1, 1, -1, cl1};

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

Line Loop(45-16) = {10, 12, -11, -9};
Line Loop(46-16) = {15, 17, -18, -14};
Plane Surface(47-16) = {46-16};
Plane Surface(48) = {45-16};
Line Loop(49) = {22, 24, -23, -21};
Plane Surface(50) = {49};
Line Loop(51) = {20, -19, -17, 18};
Plane Surface(52) = {51};
Line Loop(53) = {29, -41, -25, 44};
Plane Surface(54) = {53};
Line Loop(55) = {25, -38, -21, 39};
Plane Surface(56) = {55};
Line Loop(57) = {21, -36, -17, 33};
Plane Surface(58) = {57};
Line Loop(59) = {31, -42, -27, 41};
Plane Surface(60) = {59};
Line Loop(61) = {27, -37, -23, 38};
Plane Surface(62) = {61};
Line Loop(63) = {23, -35, -19, 36};
Plane Surface(64) = {63};
Line Loop(65) = {32, -42, -28, 43};
Plane Surface(66) = {65};
Line Loop(67) = {28, -37, -24, 40};
Plane Surface(68) = {67};
Line Loop(69) = {24, -35, -20, 34};
Plane Surface(70) = {69};
Line Loop(71) = {30, -43, -26, 44};
Plane Surface(72) = {71};
Line Loop(73) = {26, -40, -22, 39};
Plane Surface(74) = {73};
Line Loop(75) = {22, -34, -18, 33};
Plane Surface(76) = {75};
