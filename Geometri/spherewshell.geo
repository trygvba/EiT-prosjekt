R=1;
p = 5;
Rs = (1+p/100)*R;

lc1=0.18;
// Inner sphere:
//Points:
Point(1) = {0, 0, 0, lc1};
Point(2) = {R, 0, 0, lc1};
Point(3) = {-R, 0, 0, lc1};
Point(4) = {0, R, 0, lc1};
Point(5) = {0, -R, 0, lc1};
Point(6) = {0, 0, R, lc1};
Point(7) = {0, 0, -R, lc1};
//Lines:
Circle(1) = {2, 1, 4};
Circle(2) = {4,1,3};
Circle(3) = {3,1,5};
Circle(4) = {5,1,2};

Circle(5) = {6,1,5};
Circle(6) = {5,1,7};
Circle(7) = {7,1,4};
Circle(8) = {4,1,6};

Circle(9) = {2, 1, 6};
Circle(10) = {6, 1, 3};
Circle(11) = {3, 1, 7};
Circle(12) = {7, 1, 2};
//Surfaces:
Line Loop(13) = {12, -4, 6};
Ruled Surface(14) = {13};
Line Loop(15) = {12, 1, -7};
Ruled Surface(16) = {15};
Line Loop(17) = {1, 8, -9};
Ruled Surface(18) = {17};
Line Loop(19) = {9, 5, 4};
Ruled Surface(20) = {19};
Line Loop(21) = {5, -3, -10};
Ruled Surface(22) = {21};
Line Loop(23) = {3, 6, -11};
Ruled Surface(24) = {23};
Line Loop(25) = {11, 7, 2};
Ruled Surface(26) = {25};
Line Loop(27) = {10, -2, 8};
Ruled Surface(28) = {27};
//Volume:
Surface Loop(29) = {22, 20, 18, 16, 14, 24, 26, 28};
Volume(30) = {29};

//Outer Shell:
lc2=lc1;
//Points:
Point(31) = {Rs, 0, 0, lc2};
Point(32) = {-Rs, 0, 0, lc2};
Point(33) = {0, Rs, 0, lc2};
Point(34) = {0, -Rs, 0, lc2};
Point(35) = {0, 0, Rs, lc2};
Point(36) = {0, 0, -Rs, lc2};

//Lines:
Circle(31) = {31, 1, 33};
Circle(32) = {33, 1, 32};
Circle(33) = {32, 1, 34};
Circle(34) = {34, 1, 31};
Circle(35) = {34, 1, 36};
Circle(36) = {36, 1, 33};
Circle(37) = {33, 1, 35};

Circle(38) = {35, 1, 34};
Circle(39) = {32, 1, 35};
Circle(40) = {35, 1, 31};
Circle(41) = {31, 1, 36};
Circle(42) = {36, 1, 32};

//Surfaces:
Line Loop(43) = {36, 32, -42};
Ruled Surface(44) = {43};
Line Loop(45) = {36, -31, 41};
Ruled Surface(46) = {45};
Line Loop(47) = {37, 40, 31};
Ruled Surface(48) = {47};
Line Loop(49) = {32, 39, -37};
Ruled Surface(50) = {49};
Line Loop(51) = {42, 33, 35};
Ruled Surface(52) = {51};
Line Loop(53) = {35, -41, -34};
Ruled Surface(54) = {53};
Line Loop(55) = {40, -34, -38};
Ruled Surface(56) = {55};
Line Loop(57) = {39, 38, -33};
Ruled Surface(58) = {57};

//Volume:
Surface Loop(59) = {44, 46, 48, 50, 58, 56, 54, 52};
Volume(60) = {29, 59};

//Physical Groups:
Physical Volume(61) = {60}; //Outer shell
Physical Volume(62) = {30}; // Inner sphere
Physical Surface(63) = {44, 50, 48, 46, 56, 54, 52, 58}; //Outer Surface
