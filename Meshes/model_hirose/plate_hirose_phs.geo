// interplate_phs.geo

radius = 5.0;
pio2 = Pi/2;

// IMP to OK
cellsize = 0.150000;
Point(1) =  {139.313298, 35.105794, 0, cellsize};
Point(2) =  {139.332086, 35.360903, 0, cellsize};
Point(3) =  {139.362086, 35.806732, 0, cellsize};
Point(4) =  {139.367860, 36.164701, 0, cellsize};
Point(5) =  {139.046807, 36.169901, 0, cellsize};
Point(6) =  {138.816835, 36.138515, 0, cellsize};
Point(7) =  {138.618845, 36.014137, 0, cellsize};
Point(8) =  {138.471622, 35.904989, 0, cellsize};
Point(9) =  {138.314246, 35.826301, 0, cellsize};
Point(10) =  {138.404551, 35.730771, 0, cellsize};
Point(11) =  {138.441162, 35.559776, 0, cellsize};
Point(12) =  {138.451376, 35.429421, 0, cellsize};
Point(13) =  {138.487390, 35.248942, 0, cellsize};
Point(14) =  {138.642000, 35.106100, 0, cellsize};
Point(15) =  {138.680000, 35.250000, 0, cellsize};
Point(16) =  {139.000000, 35.348000, 0, cellsize};
Point(17) =  {139.188629, 35.254355, 0, cellsize};

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

// IMP to NAN
cellsize = 0.150000;
Point(18) =  {138.441162, 35.559776, 0, cellsize};
Point(19) =  {138.200000, 35.350000, 0, cellsize};
Point(20) =  {138.000000, 35.200000, 0, cellsize};
Point(21) =  {137.800000, 35.085000, 0, cellsize};
Point(22) =  {137.500000, 35.020000, 0, cellsize};
Point(23) =  {137.100000, 35.010000, 0, cellsize};
Point(24) =  {136.900000, 35.050000, 0, cellsize};
Point(25) =  {136.730000, 34.950000, 0, cellsize};
Point(26) =  {136.600000, 34.650000, 0, cellsize};
Point(27) =  {136.450000, 34.200000, 0, cellsize};
Point(28) =  {136.350000, 33.850000, 0, cellsize};
Point(29) =  {136.300000, 33.550000, 0, cellsize};
Point(30) =  {136.350000, 33.250000, 0, cellsize};
Point(31) =  {136.779000, 32.954200, 0, cellsize};
Point(32) =  {137.072224, 33.062003, 0, cellsize};
Point(33) =  {137.308701, 33.230072, 0, cellsize};
Point(34) =  {137.592330, 33.356489, 0, cellsize};
Point(35) =  {137.847472, 33.530805, 0, cellsize};
Point(36) =  {138.086035, 33.719012, 0, cellsize};
Point(37) =  {138.327724, 33.895610, 0, cellsize};
Point(38) =  {138.484017, 34.085867, 0, cellsize};
Point(39) =  {138.483129, 34.307243, 0, cellsize};
Point(40) =  {138.577802, 34.562761, 0, cellsize};
Point(41) =  {138.620987, 34.832977, 0, cellsize};
Point(42) =  {138.642000, 35.106100, 0, cellsize};
Point(43) =  {138.487390, 35.248942, 0, cellsize};
Point(44) =  {138.451376, 35.429421, 0, cellsize};

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
Line(544) = {44, 18};

Line Loop(1002) = {518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544};
Plane Surface(2002) = {1002};

// PHS to NAN
cellsize = 0.150000;
Point(45) =  {138.314246, 35.826301, 0, cellsize};
Point(46) =  {138.026525, 35.702743, 0, cellsize};
Point(47) =  {137.750101, 35.635173, 0, cellsize};
Point(48) =  {137.473677, 35.641316, 0, cellsize};
Point(49) =  {137.181896, 35.708886, 0, cellsize};
Point(50) =  {136.936186, 35.801027, 0, cellsize};
Point(51) =  {136.715046, 35.893169, 0, cellsize};
Point(52) =  {136.472408, 35.957668, 0, cellsize};
Point(53) =  {136.303482, 36.062095, 0, cellsize};
Point(54) =  {136.174484, 36.092808, 0, cellsize};
Point(55) =  {135.916488, 35.960739, 0, cellsize};
Point(56) =  {135.762919, 35.828670, 0, cellsize};
Point(57) =  {135.679992, 35.632102, 0, cellsize};
Point(58) =  {135.618564, 35.395605, 0, cellsize};
Point(59) =  {135.547922, 35.149895, 0, cellsize};
Point(60) =  {135.514137, 34.974826, 0, cellsize};
Point(61) =  {135.557137, 34.778258, 0, cellsize};
Point(62) =  {135.640064, 34.658474, 0, cellsize};
Point(63) =  {135.584779, 34.431192, 0, cellsize};
Point(64) =  {135.409710, 34.354408, 0, cellsize};
Point(65) =  {135.203928, 34.339051, 0, cellsize};
Point(66) =  {135.151715, 34.624689, 0, cellsize};
Point(67) =  {135.056502, 34.778258, 0, cellsize};
Point(68) =  {134.829220, 34.925684, 0, cellsize};
Point(69) =  {134.568153, 34.959470, 0, cellsize};
Point(70) =  {134.270229, 34.977898, 0, cellsize};
Point(71) =  {133.956948, 34.974826, 0, cellsize};
Point(72) =  {133.599000, 34.823675, 0, cellsize};
Point(73) =  {133.180032, 34.807352, 0, cellsize};
Point(74) =  {132.810035, 34.687647, 0, cellsize};
Point(75) =  {132.467244, 34.627794, 0, cellsize};
Point(76) =  {132.184305, 34.464560, 0, cellsize};
Point(77) =  {132.031953, 34.197945, 0, cellsize};
Point(78) =  {131.934012, 33.893241, 0, cellsize};
Point(79) =  {131.885042, 33.588537, 0, cellsize};
Point(80) =  {131.754455, 33.207658, 0, cellsize};
Point(81) =  {131.634750, 32.962807, 0, cellsize};
Point(82) =  {131.515045, 32.679868, 0, cellsize};
Point(83) =  {131.373575, 32.445899, 0, cellsize};
Point(84) =  {131.242988, 32.184725, 0, cellsize};
Point(85) =  {131.128724, 31.912668, 0, cellsize};
Point(86) =  {131.439206, 31.834671, 0, cellsize};
Point(87) =  {131.717290, 31.760237, 0, cellsize};
Point(88) =  {131.889634, 31.608569, 0, cellsize};
Point(89) =  {132.107530, 31.410818, 0, cellsize};
Point(90) =  {132.322425, 31.211633, 0, cellsize};
Point(91) =  {132.568319, 31.002844, 0, cellsize};
Point(92) =  {132.824000, 30.754400, 0, cellsize};
Point(93) =  {132.902806, 31.034987, 0, cellsize};
Point(94) =  {133.035121, 31.102386, 0, cellsize};
Point(95) =  {133.263259, 31.322388, 0, cellsize};
Point(96) =  {133.387000, 31.553500, 0, cellsize};
Point(97) =  {133.655840, 31.689286, 0, cellsize};
Point(98) =  {133.905723, 31.833692, 0, cellsize};
Point(99) =  {134.126687, 31.959162, 0, cellsize};
Point(100) =  {134.348578, 32.096855, 0, cellsize};
Point(101) =  {134.586050, 32.188843, 0, cellsize};
Point(102) =  {134.820578, 32.301442, 0, cellsize};
Point(103) =  {135.010961, 32.458719, 0, cellsize};
Point(104) =  {135.199000, 32.624000, 0, cellsize};
Point(105) =  {135.487550, 32.653948, 0, cellsize};
Point(106) =  {135.787922, 32.688894, 0, cellsize};
Point(107) =  {136.050270, 32.745821, 0, cellsize};
Point(108) =  {136.299633, 32.803317, 0, cellsize};
Point(109) =  {136.549827, 32.866404, 0, cellsize};
Point(110) =  {136.779000, 32.954200, 0, cellsize};
Point(111) =  {136.350000, 33.250000, 0, cellsize};
Point(112) =  {136.300000, 33.550000, 0, cellsize};
Point(113) =  {136.350000, 33.850000, 0, cellsize};
Point(114) =  {136.450000, 34.200000, 0, cellsize};
Point(115) =  {136.600000, 34.650000, 0, cellsize};
Point(116) =  {136.730000, 34.950000, 0, cellsize};
Point(117) =  {136.900000, 35.050000, 0, cellsize};
Point(118) =  {137.100000, 35.010000, 0, cellsize};
Point(119) =  {137.500000, 35.020000, 0, cellsize};
Point(120) =  {137.800000, 35.085000, 0, cellsize};
Point(121) =  {138.000000, 35.200000, 0, cellsize};
Point(122) =  {138.200000, 35.350000, 0, cellsize};
Point(123) =  {138.441162, 35.559776, 0, cellsize};
Point(124) =  {138.404551, 35.730771, 0, cellsize};

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
Line(569) = {69, 70};
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
Line(607) = {107, 108};
Line(608) = {108, 109};
Line(609) = {109, 110};
Line(610) = {110, 111};
Line(611) = {111, 112};
Line(612) = {112, 113};
Line(613) = {113, 114};
Line(614) = {114, 115};
Line(615) = {115, 116};
Line(616) = {116, 117};
Line(617) = {117, 118};
Line(618) = {118, 119};
Line(619) = {119, 120};
Line(620) = {120, 121};
Line(621) = {121, 122};
Line(622) = {122, 123};
Line(623) = {123, 124};
Line(624) = {124, 45};

Line Loop(1003) = {545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624};
Plane Surface(2003) = {1003};

// PHS to ON-North
cellsize = 0.250000;
Point(125) =  {132.824000, 30.754400, 0, cellsize};
Point(126) =  {132.568319, 31.002844, 0, cellsize};
Point(127) =  {132.322425, 31.211633, 0, cellsize};
Point(128) =  {132.107530, 31.410818, 0, cellsize};
Point(129) =  {131.889634, 31.608569, 0, cellsize};
Point(130) =  {131.717290, 31.760237, 0, cellsize};
Point(131) =  {131.439206, 31.834671, 0, cellsize};
Point(132) =  {131.226958, 31.890257, 0, cellsize};
Point(133) =  {131.118216, 31.643472, 0, cellsize};
Point(134) =  {130.966361, 31.396707, 0, cellsize};
Point(135) =  {130.781287, 31.126214, 0, cellsize};
Point(136) =  {130.600958, 30.841485, 0, cellsize};
Point(137) =  {130.449103, 30.604210, 0, cellsize};
Point(138) =  {130.243463, 30.307624, 0, cellsize};
Point(139) =  {130.041697, 30.055417, 0, cellsize};
Point(140) =  {129.848338, 29.820023, 0, cellsize};
Point(141) =  {129.646572, 29.576223, 0, cellsize};
Point(142) =  {129.436399, 29.324015, 0, cellsize};
Point(143) =  {129.234633, 29.105435, 0, cellsize};
Point(144) =  {128.969271, 28.794160, 0, cellsize};
Point(145) =  {128.688354, 28.521506, 0, cellsize};
Point(146) =  {128.498322, 28.331474, 0, cellsize};
Point(147) =  {128.250454, 28.108393, 0, cellsize};
Point(148) =  {128.440486, 27.901836, 0, cellsize};
Point(149) =  {128.605731, 27.703542, 0, cellsize};
Point(150) =  {128.812288, 27.496985, 0, cellsize};
Point(151) =  {129.010582, 27.282167, 0, cellsize};
Point(152) =  {129.233663, 27.042561, 0, cellsize};
Point(153) =  {129.431957, 26.819480, 0, cellsize};
Point(154) =  {129.669219, 26.570057, 0, cellsize};
Point(155) =  {129.960742, 26.902103, 0, cellsize};
Point(156) =  {130.167298, 27.216069, 0, cellsize};
Point(157) =  {130.365592, 27.472199, 0, cellsize};
Point(158) =  {130.537613, 27.773685, 0, cellsize};
Point(159) =  {130.764692, 27.970501, 0, cellsize};
Point(160) =  {131.008493, 28.197488, 0, cellsize};
Point(161) =  {131.253187, 28.414691, 0, cellsize};
Point(162) =  {131.479280, 28.811193, 0, cellsize};
Point(163) =  {131.709635, 29.237751, 0, cellsize};
Point(164) =  {131.966881, 29.618257, 0, cellsize};
Point(165) =  {132.232560, 30.010013, 0, cellsize};
Point(166) =  {132.546958, 30.400101, 0, cellsize};

Line(625) = {125, 126};
Line(626) = {126, 127};
Line(627) = {127, 128};
Line(628) = {128, 129};
Line(629) = {129, 130};
Line(630) = {130, 131};
Line(631) = {131, 132};
Line(632) = {132, 133};
Line(633) = {133, 134};
Line(634) = {134, 135};
Line(635) = {135, 136};
Line(636) = {136, 137};
Line(637) = {137, 138};
Line(638) = {138, 139};
Line(639) = {139, 140};
Line(640) = {140, 141};
Line(641) = {141, 142};
Line(642) = {142, 143};
Line(643) = {143, 144};
Line(644) = {144, 145};
Line(645) = {145, 146};
Line(646) = {146, 147};
Line(647) = {147, 148};
Line(648) = {148, 149};
Line(649) = {149, 150};
Line(650) = {150, 151};
Line(651) = {151, 152};
Line(652) = {152, 153};
Line(653) = {153, 154};
Line(654) = {154, 155};
Line(655) = {155, 156};
Line(656) = {156, 157};
Line(657) = {157, 158};
Line(658) = {158, 159};
Line(659) = {159, 160};
Line(660) = {160, 161};
Line(661) = {161, 162};
Line(662) = {162, 163};
Line(663) = {163, 164};
Line(664) = {164, 165};
Line(665) = {165, 166};
Line(666) = {166, 125};

Line Loop(1004) = {625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666};
Plane Surface(2004) = {1004};

// PHS to ON-Center
cellsize = 0.900000;
Point(167) =  {129.669219, 26.570057, 0, cellsize};
Point(168) =  {129.431957, 26.819480, 0, cellsize};
Point(169) =  {129.233663, 27.042561, 0, cellsize};
Point(170) =  {129.010582, 27.282167, 0, cellsize};
Point(171) =  {128.812288, 27.496985, 0, cellsize};
Point(172) =  {128.605731, 27.703542, 0, cellsize};
Point(173) =  {128.440486, 27.901836, 0, cellsize};
Point(174) =  {128.250454, 28.108393, 0, cellsize};
Point(175) =  {127.969350, 27.852334, 0, cellsize};
Point(176) =  {127.705784, 27.552414, 0, cellsize};
Point(177) =  {127.451306, 27.270671, 0, cellsize};
Point(178) =  {127.196829, 26.961663, 0, cellsize};
Point(179) =  {126.951440, 26.752628, 0, cellsize};
Point(180) =  {126.678785, 26.543592, 0, cellsize};
Point(181) =  {126.306157, 26.325469, 0, cellsize};
Point(182) =  {125.988060, 26.098257, 0, cellsize};
Point(183) =  {125.697229, 25.916487, 0, cellsize};
Point(184) =  {125.333689, 25.698363, 0, cellsize};
Point(185) =  {125.015592, 25.489328, 0, cellsize};
Point(186) =  {124.706584, 25.271205, 0, cellsize};
Point(187) =  {124.424841, 25.071258, 0, cellsize};
Point(188) =  {124.461195, 24.825869, 0, cellsize};
Point(189) =  {124.524814, 24.562303, 0, cellsize};
Point(190) =  {124.570256, 24.289648, 0, cellsize};
Point(191) =  {124.615699, 23.998817, 0, cellsize};
Point(192) =  {124.679318, 23.698897, 0, cellsize};
Point(193) =  {124.730746, 23.411937, 0, cellsize};
Point(194) =  {125.091368, 23.400519, 0, cellsize};
Point(195) =  {125.633609, 23.480773, 0, cellsize};
Point(196) =  {126.125618, 23.594746, 0, cellsize};
Point(197) =  {126.497015, 23.844312, 0, cellsize};
Point(198) =  {126.833289, 24.089701, 0, cellsize};
Point(199) =  {127.200251, 24.345569, 0, cellsize};
Point(200) =  {127.596722, 24.689541, 0, cellsize};
Point(201) =  {127.918610, 24.979507, 0, cellsize};
Point(202) =  {128.314712, 25.280293, 0, cellsize};
Point(203) =  {128.643741, 25.600289, 0, cellsize};
Point(204) =  {129.014526, 25.907399, 0, cellsize};
Point(205) =  {129.382969, 26.212710, 0, cellsize};

Line(667) = {167, 168};
Line(668) = {168, 169};
Line(669) = {169, 170};
Line(670) = {170, 171};
Line(671) = {171, 172};
Line(672) = {172, 173};
Line(673) = {173, 174};
Line(674) = {174, 175};
Line(675) = {175, 176};
Line(676) = {176, 177};
Line(677) = {177, 178};
Line(678) = {178, 179};
Line(679) = {179, 180};
Line(680) = {180, 181};
Line(681) = {181, 182};
Line(682) = {182, 183};
Line(683) = {183, 184};
Line(684) = {184, 185};
Line(685) = {185, 186};
Line(686) = {186, 187};
Line(687) = {187, 188};
Line(688) = {188, 189};
Line(689) = {189, 190};
Line(690) = {190, 191};
Line(691) = {191, 192};
Line(692) = {192, 193};
Line(693) = {193, 194};
Line(694) = {194, 195};
Line(695) = {195, 196};
Line(696) = {196, 197};
Line(697) = {197, 198};
Line(698) = {198, 199};
Line(699) = {199, 200};
Line(700) = {200, 201};
Line(701) = {201, 202};
Line(702) = {202, 203};
Line(703) = {203, 204};
Line(704) = {204, 205};
Line(705) = {205, 167};

Line Loop(1005) = {667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705};
Plane Surface(2005) = {1005};

// PHS to ON-South
cellsize = 0.900000;
Point(206) =  {124.730746, 23.411937, 0, cellsize};
Point(207) =  {124.679318, 23.698897, 0, cellsize};
Point(208) =  {124.615699, 23.998817, 0, cellsize};
Point(209) =  {124.570256, 24.289648, 0, cellsize};
Point(210) =  {124.524814, 24.562303, 0, cellsize};
Point(211) =  {124.461195, 24.825869, 0, cellsize};
Point(212) =  {124.424841, 25.071258, 0, cellsize};
Point(213) =  {124.158967, 24.918280, 0, cellsize};
Point(214) =  {123.887066, 24.815675, 0, cellsize};
Point(215) =  {123.558732, 24.743853, 0, cellsize};
Point(216) =  {123.184227, 24.718201, 0, cellsize};
Point(217) =  {122.809721, 24.713071, 0, cellsize};
Point(218) =  {122.445476, 24.713071, 0, cellsize};
Point(219) =  {122.132533, 24.718201, 0, cellsize};
Point(220) =  {121.860632, 24.713071, 0, cellsize};
Point(221) =  {121.855502, 24.543774, 0, cellsize};
Point(222) =  {121.809330, 24.374477, 0, cellsize};
Point(223) =  {121.983757, 24.148748, 0, cellsize};
Point(224) =  {122.137949, 23.948531, 0, cellsize};
Point(225) =  {122.389044, 23.815284, 0, cellsize};
Point(226) =  {122.666075, 23.661378, 0, cellsize};
Point(227) =  {122.975050, 23.479133, 0, cellsize};
Point(228) =  {123.353524, 23.456169, 0, cellsize};
Point(229) =  {123.733160, 23.445909, 0, cellsize};
Point(230) =  {124.036755, 23.437578, 0, cellsize};
Point(231) =  {124.384697, 23.420258, 0, cellsize};

Line(706) = {206, 207};
Line(707) = {207, 208};
Line(708) = {208, 209};
Line(709) = {209, 210};
Line(710) = {210, 211};
Line(711) = {211, 212};
Line(712) = {212, 213};
Line(713) = {213, 214};
Line(714) = {214, 215};
Line(715) = {215, 216};
Line(716) = {216, 217};
Line(717) = {217, 218};
Line(718) = {218, 219};
Line(719) = {219, 220};
Line(720) = {220, 221};
Line(721) = {221, 222};
Line(722) = {222, 223};
Line(723) = {223, 224};
Line(724) = {224, 225};
Line(725) = {225, 226};
Line(726) = {226, 227};
Line(727) = {227, 228};
Line(728) = {228, 229};
Line(729) = {229, 230};
Line(730) = {230, 231};
Line(731) = {231, 206};

Line Loop(1006) = {706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731};
Plane Surface(2006) = {1006};

