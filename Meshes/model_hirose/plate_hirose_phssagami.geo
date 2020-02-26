// interplate_phssagami.geo

radius = 5.0;
cellsize = 0.20;
pio2 = Pi/2;

// Sagami Trough
Point(1) =  {141.883000, 34.212700, 0, cellsize};
Point(2) =  {141.723527, 34.727051, 0, cellsize};
Point(3) =  {141.498353, 35.206267, 0, cellsize};
Point(4) =  {141.215442, 35.656616, 0, cellsize};
Point(5) =  {140.753546, 35.951074, 0, cellsize};
Point(6) =  {140.418672, 36.297495, 0, cellsize};
Point(7) =  {140.291651, 36.436064, 0, cellsize};
Point(8) =  {140.043382, 36.343685, 0, cellsize};
Point(9) =  {139.754697, 36.280174, 0, cellsize};
Point(10) =  {139.367860, 36.164701, 0, cellsize};
Point(11) =  {139.362086, 35.806732, 0, cellsize};
Point(12) =  {139.332086, 35.360903, 0, cellsize};
Point(13) =  {139.338000, 34.920000, 0, cellsize};
Point(14) =  {139.872122, 34.727889, 0, cellsize};
Point(15) =  {140.403744, 34.637873, 0, cellsize};
Point(16) =  {140.901336, 34.517400, 0, cellsize};
Point(17) =  {141.407766, 34.364671, 0, cellsize};

Line(501) = {1, 2};
Line(502) = {2, 3};
Line(503) = {3, 4};
Line(504) = {4, 5};
Line(505) = {5, 6};
Line(506) = {6, 7};
Line(507) = {7, 8};
Line(508) = {8, 9};
Line(509) = {9, 10};
Line(510) = {10, 11};
Line(511) = {11, 12};
Line(512) = {12, 13};
Line(513) = {13, 14};
Line(514) = {14, 15};
Line(515) = {15, 16};
Line(516) = {16, 17};
Line(517) = {17, 1};

Line Loop(1001) = {501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517};
Plane Surface(2001) = {1001};

