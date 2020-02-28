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

// IMP to THE
Point(20) =  {139.332086, 35.360903, 0, cellsize};
Point(21) =  {139.364963, 35.858868, 0, cellsize};
Point(22) =  {139.364963, 36.291443, 0, cellsize};
Point(23) =  {139.186652, 36.298420, 0, cellsize};
Point(24) =  {139.043081, 36.139466, 0, cellsize};
Point(25) =  {139.058464, 35.862579, 0, cellsize};
Point(26) =  {139.027699, 35.570309, 0, cellsize};
Point(27) =  {139.000000, 35.348000, 0, cellsize};
Point(28) =  {139.188629, 35.254355, 0, cellsize};
Point(29) =  {139.313298, 35.105794, 0, cellsize};

Line(520) = {20, 21};
Line(521) = {21, 22};
Line(522) = {22, 23};
Line(523) = {23, 24};
Line(524) = {24, 25};
Line(525) = {25, 26};
Line(526) = {26, 27};
Line(527) = {27, 28};
Line(528) = {28, 29};
Line(529) = {29, 20};

Line Loop(1002) = {520, 521, 522, 523, 524, 525, 526, 527, 528, 529};
Plane Surface(2002) = {1002};

// IMP to THW
Point(30) =  {138.642000, 35.106100, 0, cellsize};
Point(31) =  {138.680000, 35.250000, 0, cellsize};
Point(32) =  {139.000000, 35.348000, 0, cellsize};
Point(33) =  {139.027699, 35.570309, 0, cellsize};
Point(34) =  {139.058464, 35.862579, 0, cellsize};
Point(35) =  {139.043081, 36.139466, 0, cellsize};
Point(36) =  {139.186652, 36.298420, 0, cellsize};
Point(37) =  {138.840845, 36.268802, 0, cellsize};
Point(38) =  {138.543304, 36.195721, 0, cellsize};
Point(39) =  {138.329282, 36.096541, 0, cellsize};
Point(40) =  {138.143341, 36.001782, 0, cellsize};
Point(41) =  {138.275481, 35.866585, 0, cellsize};
Point(42) =  {138.404551, 35.730771, 0, cellsize};
Point(43) =  {138.437497, 35.587068, 0, cellsize};
Point(44) =  {138.451376, 35.429421, 0, cellsize};
Point(45) =  {138.487390, 35.248942, 0, cellsize};

Line(530) = {30, 31};
Line(531) = {31, 32};
Line(532) = {32, 33};
Line(533) = {33, 34};
Line(534) = {34, 35};
Line(535) = {35, 36};
Line(536) = {36, 37};
Line(537) = {37, 38};
Line(538) = {38, 39};
Line(539) = {39, 40};
Line(540) = {40, 41};
Line(541) = {41, 42};
Line(542) = {42, 43};
Line(543) = {43, 44};
Line(544) = {44, 45};
Line(545) = {45, 30};

Line Loop(1003) = {530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545};
Plane Surface(2003) = {1003};

