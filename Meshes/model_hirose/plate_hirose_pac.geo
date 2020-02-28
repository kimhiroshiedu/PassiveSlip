// interplate_pac.geo

radius = 5.0;
cellsize = 0.15;
pio2 = Pi/2;

// Kuril trench
Point(1) =  {144.982390, 41.171590, 0, cellsize};
Point(2) =  {145.569683, 41.602876, 0, cellsize};
Point(3) =  {146.115016, 42.027286, 0, cellsize};
Point(4) =  {146.779507, 42.400754, 0, cellsize};
Point(5) =  {147.449021, 42.745782, 0, cellsize};
Point(6) =  {148.095585, 43.084267, 0, cellsize};
Point(7) =  {148.729320, 43.381928, 0, cellsize};
Point(8) =  {149.387088, 43.737738, 0, cellsize};
Point(9) =  {150.115045, 44.080193, 0, cellsize};
Point(10) =  {150.749607, 44.395449, 0, cellsize};
Point(11) =  {150.455961, 44.902209, 0, cellsize};
Point(12) =  {150.150569, 45.408970, 0, cellsize};
Point(13) =  {149.944794, 45.877731, 0, cellsize};
Point(14) =  {149.501433, 45.621945, 0, cellsize};
Point(15) =  {148.938705, 45.315003, 0, cellsize};
Point(16) =  {148.495343, 44.854589, 0, cellsize};
Point(17) =  {148.051982, 44.445332, 0, cellsize};
Point(18) =  {147.421045, 44.206599, 0, cellsize};
Point(19) =  {146.738950, 44.070180, 0, cellsize};
Point(20) =  {146.125065, 43.831447, 0, cellsize};
Point(21) =  {145.528232, 43.558609, 0, cellsize};
Point(22) =  {144.931400, 43.302824, 0, cellsize};
Point(23) =  {144.334567, 43.098196, 0, cellsize};
Point(24) =  {143.703630, 42.944724, 0, cellsize};
Point(25) =  {143.072692, 42.910620, 0, cellsize};
Point(26) =  {142.697540, 42.808305, 0, cellsize};
Point(27) =  {143.036685, 42.316727, 0, cellsize};
Point(28) =  {143.443646, 41.738414, 0, cellsize};
Point(29) =  {144.043379, 41.224357, 0, cellsize};
Point(30) =  {144.461120, 40.773660, 0, cellsize};

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
Line(519) = {19, 20};
Line(520) = {20, 21};
Line(521) = {21, 22};
Line(522) = {22, 23};
Line(523) = {23, 24};
Line(524) = {24, 25};
Line(525) = {25, 26};
Line(526) = {26, 27};
Line(527) = {27, 28};
Line(528) = {28, 29};
Line(529) = {29, 30};
Line(530) = {30, 1};

Line Loop(1001) = {501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530};
Plane Surface(2001) = {1001};

// Japan trench
Point(31) =  {141.883000, 34.212700, 0, cellsize};
Point(32) =  {142.120000, 34.830000, 0, cellsize};
Point(33) =  {142.340000, 35.400000, 0, cellsize};
Point(34) =  {142.539427, 35.926512, 0, cellsize};
Point(35) =  {142.870000, 36.320000, 0, cellsize};
Point(36) =  {143.290000, 36.860000, 0, cellsize};
Point(37) =  {143.609750, 37.250333, 0, cellsize};
Point(38) =  {143.810000, 37.870000, 0, cellsize};
Point(39) =  {143.990000, 38.460000, 0, cellsize};
Point(40) =  {144.140604, 38.981194, 0, cellsize};
Point(41) =  {144.250000, 39.580000, 0, cellsize};
Point(42) =  {144.360000, 40.200000, 0, cellsize};
Point(43) =  {144.461120, 40.773660, 0, cellsize};
Point(44) =  {144.043379, 41.224357, 0, cellsize};
Point(45) =  {143.443646, 41.738414, 0, cellsize};
Point(46) =  {143.036685, 42.316727, 0, cellsize};
Point(47) =  {142.697540, 42.808305, 0, cellsize};
Point(48) =  {142.237131, 42.484310, 0, cellsize};
Point(49) =  {141.810822, 42.211473, 0, cellsize};
Point(50) =  {141.520932, 42.023897, 0, cellsize};
Point(51) =  {141.401566, 41.478221, 0, cellsize};
Point(52) =  {141.374612, 40.937184, 0, cellsize};
Point(53) =  {141.310387, 40.335073, 0, cellsize};
Point(54) =  {141.238133, 39.684794, 0, cellsize};
Point(55) =  {141.141796, 39.138881, 0, cellsize};
Point(56) =  {141.037430, 38.657192, 0, cellsize};
Point(57) =  {140.892923, 38.199588, 0, cellsize};
Point(58) =  {140.676164, 37.653675, 0, cellsize};
Point(59) =  {140.411236, 37.155930, 0, cellsize};
Point(60) =  {140.023846, 36.573932, 0, cellsize};
Point(61) =  {140.286737, 36.433882, 0, cellsize};
Point(62) =  {140.485896, 36.211451, 0, cellsize};
Point(63) =  {140.740819, 35.956067, 0, cellsize};
Point(64) =  {140.947945, 35.824256, 0, cellsize};
Point(65) =  {141.202869, 35.667730, 0, cellsize};
Point(66) =  {141.354230, 35.437060, 0, cellsize};
Point(67) =  {141.521523, 35.173438, 0, cellsize};
Point(68) =  {141.672884, 34.852148, 0, cellsize};
Point(69) =  {141.784413, 34.563811, 0, cellsize};

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
Line(545) = {45, 46};
Line(546) = {46, 47};
Line(547) = {47, 48};
Line(548) = {48, 49};
Line(549) = {49, 50};
Line(550) = {50, 51};
Line(551) = {51, 52};
Line(552) = {52, 53};
Line(553) = {53, 54};
Line(554) = {54, 55};
Line(555) = {55, 56};
Line(556) = {56, 57};
Line(557) = {57, 58};
Line(558) = {58, 59};
Line(559) = {59, 60};
Line(560) = {60, 61};
Line(561) = {61, 62};
Line(562) = {62, 63};
Line(563) = {63, 64};
Line(564) = {64, 65};
Line(565) = {65, 66};
Line(566) = {66, 67};
Line(567) = {67, 68};
Line(568) = {68, 69};
Line(569) = {69, 31};

Line Loop(1002) = {531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569};
Plane Surface(2002) = {1002};

// Izu-Ogasawara trench
Point(70) =  {142.288514, 30.545359, 0, cellsize};
Point(71) =  {142.235935, 30.942721, 0, cellsize};
Point(72) =  {142.179777, 31.347760, 0, cellsize};
Point(73) =  {142.130000, 31.730000, 0, cellsize};
Point(74) =  {142.091528, 32.054737, 0, cellsize};
Point(75) =  {142.034438, 32.337184, 0, cellsize};
Point(76) =  {142.019325, 32.702799, 0, cellsize};
Point(77) =  {142.006957, 33.059370, 0, cellsize};
Point(78) =  {141.987235, 33.417141, 0, cellsize};
Point(79) =  {141.931077, 33.807451, 0, cellsize};
Point(80) =  {141.883000, 34.212700, 0, cellsize};
Point(81) =  {141.784413, 34.563811, 0, cellsize};
Point(82) =  {141.672884, 34.852148, 0, cellsize};
Point(83) =  {141.521523, 35.173438, 0, cellsize};
Point(84) =  {141.354230, 35.437060, 0, cellsize};
Point(85) =  {141.202869, 35.667730, 0, cellsize};
Point(86) =  {140.947945, 35.824256, 0, cellsize};
Point(87) =  {140.740819, 35.956067, 0, cellsize};
Point(88) =  {140.485896, 36.211451, 0, cellsize};
Point(89) =  {140.286737, 36.433882, 0, cellsize};
Point(90) =  {140.023846, 36.573932, 0, cellsize};
Point(91) =  {139.833210, 36.176497, 0, cellsize};
Point(92) =  {139.801098, 35.710865, 0, cellsize};
Point(93) =  {139.785041, 35.205092, 0, cellsize};
Point(94) =  {139.872122, 34.727889, 0, cellsize};
Point(95) =  {139.985143, 34.248304, 0, cellsize};
Point(96) =  {140.119717, 33.836664, 0, cellsize};
Point(97) =  {140.183047, 33.432940, 0, cellsize};
Point(98) =  {140.285957, 33.029217, 0, cellsize};
Point(99) =  {140.404699, 32.680906, 0, cellsize};
Point(100) =  {140.499693, 32.229685, 0, cellsize};
Point(101) =  {140.578854, 31.778465, 0, cellsize};
Point(102) =  {140.642183, 31.382657, 0, cellsize};
Point(103) =  {140.650099, 31.018514, 0, cellsize};
Point(104) =  {140.760926, 30.670204, 0, cellsize};
Point(105) =  {140.871752, 30.424803, 0, cellsize};
Point(106) =  {141.338805, 30.472300, 0, cellsize};
Point(107) =  {141.908768, 30.519797, 0, cellsize};

Line(570) = {70, 71};
Line(571) = {71, 72};
Line(572) = {72, 73};
Line(573) = {73, 74};
Line(574) = {74, 75};
Line(575) = {75, 76};
Line(576) = {76, 77};
Line(577) = {77, 78};
Line(578) = {78, 79};
Line(579) = {79, 80};
Line(580) = {80, 81};
Line(581) = {81, 82};
Line(582) = {82, 83};
Line(583) = {83, 84};
Line(584) = {84, 85};
Line(585) = {85, 86};
Line(586) = {86, 87};
Line(587) = {87, 88};
Line(588) = {88, 89};
Line(589) = {89, 90};
Line(590) = {90, 91};
Line(591) = {91, 92};
Line(592) = {92, 93};
Line(593) = {93, 94};
Line(594) = {94, 95};
Line(595) = {95, 96};
Line(596) = {96, 97};
Line(597) = {97, 98};
Line(598) = {98, 99};
Line(599) = {99, 100};
Line(600) = {100, 101};
Line(601) = {101, 102};
Line(602) = {102, 103};
Line(603) = {103, 104};
Line(604) = {104, 105};
Line(605) = {105, 106};
Line(606) = {106, 107};
Line(607) = {107, 70};

Line Loop(1003) = {570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607};
Plane Surface(2003) = {1003};

