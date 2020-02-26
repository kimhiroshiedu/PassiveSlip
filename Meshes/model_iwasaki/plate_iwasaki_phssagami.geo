// interplate_phssagami.geo

radius = 5.0;
cellsize = 0.15;
pio2 = Pi/2;

// Sagami Trough
Point(1) =  {139.338000, 34.920000, 0, cellsize};
Point(2) =  {139.872122, 34.727889, 0, cellsize};
Point(3) =  {140.403744, 34.637873, 0, cellsize};
Point(4) =  {140.901336, 34.517400, 0, cellsize};
Point(5) =  {141.407766, 34.364671, 0, cellsize};
Point(6) =  {141.883000, 34.212700, 0, cellsize};
Point(7) =  {141.677325, 34.757310, 0, cellsize};
Point(8) =  {141.370471, 35.345813, 0, cellsize};
Point(9) =  {140.932109, 36.135515, 0, cellsize};
Point(10) =  {140.784162, 36.371922, 0, cellsize};
Point(11) =  {140.685530, 36.593240, 0, cellsize};
Point(12) =  {140.581419, 36.769288, 0, cellsize};
Point(13) =  {140.438951, 36.593240, 0, cellsize};
Point(14) =  {140.219770, 36.432282, 0, cellsize};
Point(15) =  {139.945794, 36.326653, 0, cellsize};
Point(16) =  {139.633460, 36.291443, 0, cellsize};
Point(17) =  {139.364963, 36.291443, 0, cellsize};
Point(18) =  {139.364963, 35.858868, 0, cellsize};
Point(19) =  {139.332086, 35.360903, 0, cellsize};

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
Line(517) = {17, 18};
Line(518) = {18, 19};
Line(519) = {19, 1};

Line Loop(1001) = {501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519};
Plane Surface(2001) = {1001};
