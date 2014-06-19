L1 = 5;
L2 = 2;

cl1 = 0.4;
cl2 = 0.3;
Point(1) = {-1, -1, L1, cl1};
Point(2) = {-1, 1, L1, cl1};
Point(3) = {1, -1, L1, cl1};
Point(4) = {1, 1, L1, cl1};
Point(5) = {-1, -1, L2, cl2};
Point(6) = {-1, 1, L2, cl2};
Point(7) = {1, -1, L2, cl2};
Point(8) = {1, 1, L2, cl2};
Point(9) = {-1, -1, -L2, cl2};
Point(10) = {-1, 1, -L2, cl2};
Point(11) = {1, -1, -L2, cl2};
Point(12) = {1, 1, -L2, cl2};
Point(13) = {-1, -1, -L1, cl1};
Point(14) = {-1, 1, -L1, cl1};
Point(15) = {1, -1, -L1, cl1};
Point(16) = {1, 1, -L1, cl1};
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


// SPHERE PART

R=0.9;
p = 0.1;
Rs = (1+p)*R;

lc1=0.25;
// Inner sphere:
//Points:
Point(73) = {0, 0, 0, lc1};
Point(74) = {R, 0, 0, lc1};
Point(75) = {-R, 0, 0, lc1};
Point(76) = {0, R, 0, lc1};
Point(77) = {0, -R, 0, lc1};
Point(78) = {0, 0, R, lc1};
Point(79) = {0, 0, -R, lc1};
//Lines:
Circle(73) = {74, 73, 76};
Circle(74) = {76,73,75};
Circle(75) = {75,73,77};
Circle(76) = {77,73,74};

Circle(77) = {78,73,77};
Circle(78) = {77,73,79};
Circle(79) = {79,73,76};
Circle(80) = {76,73,78};

Circle(81) = {74, 73, 78};
Circle(82) = {78, 73, 75};
Circle(83) = {75, 73, 79};
Circle(84) = {79, 73, 74};
//Surfaces:
Line Loop(85) = {84, -76, 78};
Ruled Surface(86) = {85};
Line Loop(87) = {84, 73, -79};
Ruled Surface(88) = {87};
Line Loop(89) = {73, 80, -81};
Ruled Surface(90) = {89};
Line Loop(91) = {81, 77, 76};
Ruled Surface(92) = {91};
Line Loop(93) = {77, -75, -82};
Ruled Surface(94) = {93};
Line Loop(95) = {75, 78, -83};
Ruled Surface(96) = {95};
Line Loop(97) = {83, 79, 74};
Ruled Surface(98) = {97};
Line Loop(99) = {82, -74, 80};
Ruled Surface(100) = {99};
//Volume:
Surface Loop(101) = {94, 92, 90, 88, 86, 96, 98, 100};
Volume(102) = {101};

//Outer Shell:
lc2=0.2;
//Points:
Point(103) = {Rs, 0, 0, lc2};
Point(104) = {-Rs, 0, 0, lc2};
Point(105) = {0, Rs, 0, lc2};
Point(106) = {0, -Rs, 0, lc2};
Point(107) = {0, 0, Rs, lc2};
Point(108) = {0, 0, -Rs, lc2};

//Lines:
Circle(103) = {103, 73, 105};
Circle(104) = {105, 73, 104};
Circle(105) = {104, 73, 106};
Circle(106) = {106, 73, 103};
Circle(107) = {106, 73, 108};
Circle(108) = {108, 73, 105};
Circle(109) = {105, 73, 107};

Circle(110) = {107, 73, 106};
Circle(111) = {104, 73, 107};
Circle(112) = {107, 73, 103};
Circle(113) = {103, 73, 108};
Circle(114) = {108, 73, 104};

//Surfaces:
Line Loop(115) = {108, 104, -114};
Ruled Surface(116) = {115};
Line Loop(117) = {108, -103, 113};
Ruled Surface(118) = {117};
Line Loop(119) = {109, 112, 103};
Ruled Surface(120) = {119};
Line Loop(121) = {104, 111, -109};
Ruled Surface(122) = {121};
Line Loop(123) = {114, 105, 107};
Ruled Surface(124) = {123};
Line Loop(125) = {107, -113, -106};
Ruled Surface(126) = {125};
Line Loop(127) = {112, -106, -110};
Ruled Surface(128) = {127};
Line Loop(129) = {111, 110, -105};
Ruled Surface(130) = {129};

//Volume:
Surface Loop(131) = {116, 118, 120, 122, 130, 128, 126, 124};
Volume(132) = {101, 131};

//Physical Groups:
//Spheres

Physical Volume(133) = {102}; // Inner sphere
Physical Volume(134) = {132}; // Outer shell


// Physical Surface(135) = {116, 122, 120, 118, 128, 126, 124, 130}; //Outer Surface

Surface Loop(63) = {44, 48, 50, 52, 46, 54}; // Top region
Volume(64) = {63};
Surface Loop(65) = {32, 54, 58, 60, 62, 56}; // Middle region
Volume(66) = {131, 65};
Surface Loop(67) = {34, 38, 30, 40, 36, 32}; // Bottom region
Volume(68) = {67};

Physical Surface(69) = {40, 36, 30, 38, 34, 58, 60, 62, 56, 50, 48, 44, 52, 46}; // Full boundary

Physical Volume(70) = {64}; // Region 1, top region
Physical Volume(71) = {66}; // Region 2, middle region
Physical Volume(72) = {68}; // Region 3, bottom region


