/*
static const unsigned int NB = 285;//15 + 1 + 11 + 1 + 18 + 1 + 58 + 1 + 34 + 1 + 42 + 1 + 49 + 1 + 12 + 1 + 21 + 1 + 2 + 1 + 12 + 1;
static const unsigned int ls_bound[NB]=
{117, 161, 8, 48, 38, 80, 122, 166, 211, 18, 68, 118, 46, 100, 154, 210, 1, 35, 76, 120, 162, 206, 251, 296, 342, 389, 436, 485, 82, 28, 63, 99, 134, 173, 210, 37, 76, 5, 47, 90, 138, 183, 228, 274, 321, 370, 420, 76, 108, 146, 180, 213, 247, 282, 317, 353, 389, 426, 464, 502, 542, 585, 626, 668, 711, 754, 799, 845, 891, 938, 987, 1036, 1086, 1137, 1190, 1243, 1297, 1352, 1408, 1465, 1527, 1588, 1651, 1715, 1781, 1847, 10, 81, 153, 226, 300, 377, 455, 534, 616, 699, 784, 870, 959, 1053, 1148, 1246, 1346, 1448, 1553, 1660, 79, 32, 67, 102, 137, 174, 212, 249, 288, 41, 81, 122, 167, 211, 256, 304, 351, 398, 446, 40, 91, 143, 196, 249, 303, 359, 416, 474, 533, 596, 659, 721, 785, 854, 924, 80, 111, 144, 181, 215, 250, 284, 319, 355, 392, 429, 467, 507, 546, 586, 627, 670, 713, 758, 803, 27, 74, 122, 171, 222, 275, 328, 380, 434, 489, 547, 605, 664, 724, 786, 848, 912, 18, 85, 154, 225, 296, 370, 82, 113, 39, 71, 106, 142, 176, 210, 245, 281, 318, 355, 392, 430, 469, 509, 550, 591, 633, 676, 720, 764, 809, 855, 902, 951, 1000, 1050, 1100, 1151, 1204, 1258, 1314, 1370, 1428, 1487, 1547, 1608, 1670, 1733, 1798, 1864, 1932, 2001, 2072, 2234, 2312, 2391, 105, 196, 83, 123, 157, 191, 224, 258, 292, 330, 366, 403, 440, 478, 517, 76, 109, 150, 184, 218, 254, 290, 347, 385, 424, 464, 504, 545, 587, 630, 673, 717, 762, 808, 855, 904, 955, 145, 197, 250, 1, 41, 3, 30, 74, 119, 165, 212, 260, 309, 359, 32, 277};
static const unsigned int run_bound[NB]=
{282917, 282917, 282918, 282918, 282919, 282919, 282919, 282919, 282919, 282922, 282922, 282922, 282923, 282923, 282923, 282923, 283042, 283043, 283043, 283043, 283043, 283043, 283043, 283043, 283043, 283043, 283043, 283043, 283049, 283050, 283050, 283050, 283050, 283050, 283050, 283052, 283052, 283059, 283059, 283059, 283059, 283059, 283059, 283059, 283059, 283059, 283059, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283270, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283283, 283305, 283306, 283306, 283306, 283306, 283306, 283306, 283306, 283306, 283307, 283307, 283307, 283307, 283307, 283307, 283307, 283307, 283307, 283307, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283308, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283353, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283358, 283359, 283359, 283359, 283359, 283359, 283359, 283407, 283407, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283408, 283416, 283416, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283453, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283478, 283548, 283548, 283548, 283680, 283680, 283681, 283682, 283682, 283682, 283682, 283682, 283682, 283682, 283682, 283685, 283685};
*/

static const unsigned int NB = 1035;//20 + 1 + 9 + 1 + 3 + 1 + 11 + 1 + 41 + 1 + 3 + 1 + 8 + 1 + 2 + 1 + 1 + 1 + 17 + 1 + 23 + 1 + 38 + 1 + 20 + 1 + 7 + 1 + 41 + 1 + 38 + 1 + 51 + 1 + 22 + 1 + 8 + 1 + 37 + 1 + 5 + 1 + 7 + 1 + 48 + 1 + 9 + 1 + 24 + 1 + 49 + 1 + 29 + 1 + 44 + 1 + 51 + 1 + 37 + 1 + 1 + 1 + 5 + 1 + 38 + 1 + 47 + 1 + 40 + 1 + 36 + 1 + 30 + 1 + 1 + 1 + 47 + 1 + 44 + 1 + 2 + 1;
static const unsigned int ls_bound[NB]=
{70, 45, 99, 153, 209, 266, 328, 387, 448, 33, 96, 161, 229, 296, 364, 434, 505, 579, 657, 736, 82, 55, 99, 138, 178, 219, 259, 300, 342, 384, 428, 79, 117, 162, 202, 68, 121, 162, 204, 247, 291, 335, 379, 17, 64, 110, 157, 70, 115, 154, 193, 233, 274, 316, 358, 401, 444, 488, 536, 582, 628, 676, 724, 773, 822, 872, 925, 978, 1030, 1085, 1139, 1194, 1250, 1309, 1367, 1428, 27, 90, 154, 220, 286, 353, 421, 498, 569, 641, 714, 791, 18, 71, 116, 154, 193, 82, 140, 199, 259, 319, 380, 62, 132, 199, 71, 147, 221, 68, 182, 118, 185, 241, 289, 333, 375, 417, 460, 503, 61, 107, 154, 202, 251, 300, 350, 401, 453, 100, 140, 186, 228, 271, 315, 359, 405, 451, 498, 546, 601, 650, 700, 751, 803, 856, 910, 964, 1019, 1078, 1137, 1196, 1256, 77, 120, 160, 200, 247, 27, 69, 113, 157, 203, 248, 294, 341, 389, 439, 488, 538, 589, 641, 693, 746, 800, 855, 911, 968, 1026, 1084, 1143, 1205, 30, 55, 125, 196, 269, 343, 419, 495, 573, 653, 68, 113, 152, 193, 234, 275, 317, 360, 404, 448, 493, 538, 584, 632, 680, 729, 778, 829, 880, 932, 986, 77, 14, 23, 17, 59, 102, 146, 190, 71, 2, 40, 79, 119, 160, 201, 242, 284, 327, 371, 415, 460, 506, 554, 602, 651, 700, 750, 802, 855, 909, 963, 1018, 1074, 1131, 1189, 1248, 1308, 1370, 1434, 1499, 1565, 1632, 1700, 1769, 1839, 1911, 1984, 2059, 2135, 2215, 71, 116, 156, 197, 240, 284, 329, 372, 416, 460, 550, 599, 647, 6, 56, 107, 159, 212, 265, 319, 374, 430, 487, 546, 605, 666, 726, 789, 853, 933, 999, 1067, 1135, 1206, 1277, 1350, 1424, 1500, 1577, 68, 111, 148, 187, 226, 266, 306, 347, 388, 431, 474, 518, 563, 615, 661, 708, 28, 78, 129, 180, 231, 285, 338, 392, 447, 504, 562, 621, 681, 742, 804, 867, 930, 994, 1059, 1125, 1193, 1261, 1330, 1407, 1480, 1554, 1630, 25, 104, 184, 266, 348, 432, 517, 604, 693, 77, 119, 157, 196, 236, 277, 318, 359, 401, 444, 488, 533, 578, 624, 671, 719, 767, 816, 866, 917, 968, 1020, 1072, 61, 106, 144, 183, 223, 264, 304, 345, 387, 75, 119, 169, 214, 258, 302, 347, 392, 445, 491, 537, 582, 627, 672, 718, 767, 815, 864, 919, 971, 1023, 1076, 1130, 1188, 1245, 1303, 1362, 1422, 1485, 1547, 1610, 1674, 1741, 1807, 1875, 1943, 2012, 2084, 72, 116, 153, 191, 229, 268, 79, 124, 164, 205, 250, 293, 337, 381, 84, 123, 163, 210, 251, 293, 336, 380, 424, 468, 513, 559, 606, 656, 705, 756, 807, 858, 910, 963, 1016, 1071, 1127, 1183, 1240, 1299, 1359, 1420, 1481, 1546, 1612, 1678, 1746, 1815, 1884, 1954, 2027, 2100, 2174, 2249, 2325, 2403, 2485, 2567, 2650, 2736, 2822, 2910, 2999, 79, 115, 156, 192, 230, 268, 306, 346, 386, 427, 70, 106, 142, 183, 220, 259, 298, 336, 374, 413, 453, 493, 534, 576, 619, 663, 708, 777, 823, 870, 918, 966, 1015, 1065, 1116, 85, 121, 7, 44, 16, 63, 108, 154, 207, 255, 303, 352, 402, 453, 505, 562, 37, 92, 148, 57, 114, 172, 232, 292, 354, 416, 480, 544, 1, 67, 134, 204, 278, 352, 428, 504, 581, 665, 746, 828, 911, 995, 1080, 1168, 1256, 24, 120, 217, 316, 422, 4, 41, 25, 63, 101, 140, 180, 221, 19, 62, 106, 151, 196, 242, 289, 337, 385, 435, 485, 536, 588, 640, 694, 749, 805, 
861, 85, 35, 152, 221, 1, 38, 75, 113, 153, 193, 234, 276, 318, 360, 403, 450, 494, 540, 586, 7, 55, 104, 154, 205, 256, 308, 360, 414, 468, 525, 581, 638, 697, 756, 817, 878, 943, 1006, 1072, 1138, 1205, 1274, 1345, 1417, 58, 132, 209, 286, 369, 49, 5, 44, 84, 124, 165, 207, 250, 293, 337, 384, 429, 475, 522, 570, 618, 667, 718, 770, 822, 19, 72, 126, 181, 238, 309, 368, 428, 504, 568, 268, 51, 118, 186, 262, 334, 49, 122, 204, 280, 358, 438, 519, 602, 689, 774, 861, 949, 1039, 1148, 1242, 1340, 64, 7, 47, 87, 127, 168, 210, 254, 298, 342, 387, 433, 480, 530, 583, 634, 695, 746, 797, 850, 913, 968, 1023, 1080, 1138, 1199, 1258, 1318, 1379, 1442, 1505, 1572, 1637, 1713, 1782, 1852, 1922, 1993, 101, 900, 1, 71, 6, 80, 156, 237, 1, 40, 80, 121, 162, 204, 247, 291, 335, 381, 427, 475, 524, 573, 623, 674, 727, 780, 834, 889, 945, 1002, 1060, 1119, 1179, 1242, 1304, 1368, 1433, 1499, 1567, 1636, 1706, 1779, 1854, 1930, 2007, 2091, 2172, 99, 140, 175, 212, 249, 286, 325, 365, 405, 446, 488, 532, 576, 621, 666, 713, 761, 811, 860, 911, 962, 1031, 20, 75, 132, 189, 247, 20, 82, 144, 209, 274, 341, 411, 481, 553, 627, 702, 780, 859, 940, 1023, 1108, 1196, 1285, 1376, 1469, 1570, 125, 159, 193, 228, 269, 305, 343, 380, 420, 459, 499, 540, 582, 625, 670, 714, 759, 806, 853, 902, 950, 1001, 1053, 1106, 1159, 1214, 1270, 1327, 1384, 1441, 1500, 1560, 1621, 1685, 1749, 1813, 1879, 1947, 2016, 2086, 2157, 82, 6, 11, 14, 55, 97, 137, 178, 220, 263, 307, 352, 397, 444, 492, 541, 590, 640, 692, 745, 805, 859, 915, 972, 1031, 1091, 1152, 1214, 1277, 1342, 1410, 1478, 1548, 1619, 1691, 1765, 1841, 92, 127, 164, 205, 651, 698, 745, 794, 845, 896, 948, 1001, 1056, 1112, 1168, 1226, 1286, 1346, 1407, 1469, 1533, 1598, 1664, 1732, 1801, 1872, 1945, 2021, 2097, 2175, 2260, 1, 35, 89, 127, 163, 33, 72, 107, 143, 8, 45, 15, 54, 94, 136, 1, 44, 88, 133, 179, 225, 273, 322, 44, 94, 146, 199, 253, 308, 364, 422, 480, 539, 599, 664, 726, 790, 855, 921, 989, 1058, 1129, 1202, 1280, 1357, 1435, 1516, 1598, 1683, 1774, 1, 40, 76, 113, 151, 191, 232, 275, 318, 358, 23, 66, 110, 154, 199, 245, 293, 22, 70, 119, 170, 222, 275, 329, 386, 443, 502, 561, 621, 682, 744, 808, 873, 939, 1006, 1075, 1147, 1220, 1294, 1370, 1448, 1528, 1614, 1697, 1782, 1, 33, 67};
static const unsigned int run_bound[NB]=
{278873, 278874, 278874, 278874, 278874, 278874, 278874, 278874, 278874, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278875, 278923, 278923, 278923, 278923, 278923, 278923, 278923, 278923, 278923, 278923, 278957, 278957, 278957, 278957, 278962, 278962, 278962, 278962, 278962, 278962, 278962, 278962, 278963, 278963, 278963, 278963, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278969, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278975, 278976, 278986, 278986, 278986, 278986, 279024, 279024, 279024, 279024, 279024, 279024, 279029, 279029, 279029, 279071, 279071, 279071, 279080, 279080, 279115, 279115, 279115, 279115, 279115, 279115, 279115, 279115, 279115, 279116, 279116, 279116, 279116, 279116, 279116, 279116, 279116, 279116, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279588, 279653, 279653, 279653, 279653, 279653, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279654, 279656, 279658, 279658, 279658, 279658, 279658, 279658, 279658, 279658, 279658, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279667, 279681, 279682, 279683, 279685, 279685, 279685, 279685, 279685, 279691, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279694, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279715, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279716, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279760, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279766, 279767, 279767, 279767, 279767, 279767, 279767, 279767, 279767, 279767, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279794, 279823, 279823, 
279823, 279823, 279823, 279823, 279823, 279823, 279823, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279841, 279844, 279844, 279844, 279844, 279844, 279844, 279887, 279887, 279887, 279887, 279887, 279887, 279887, 279887, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279931, 279966, 279966, 279966, 279966, 279966, 279966, 279966, 279966, 279966, 279966, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279975, 279993, 279993, 279994, 279994, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280015, 280016, 280016, 280016, 280017, 280017, 280017, 280017, 280017, 280017, 280017, 280017, 280017, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280018, 280024, 280024, 280024, 280024, 280024, 280187, 280187, 280188, 280188, 280188, 280188, 280188, 280188, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280191, 280194, 280194, 280194, 280194, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280242, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280249, 280251, 280251, 280251, 280251, 280251, 280327, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280330, 280349, 280349, 280349, 280349, 280349, 280349, 280349, 280349, 280349, 280349, 280349, 280363, 280363, 280363, 280363, 280363, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280364, 280383, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 280385, 281613, 281613, 281639, 281639, 281641, 281641, 281641, 281641, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693,
281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281693, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281707, 281726, 281726, 281726, 281726, 281726, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281727, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 281797, 282033, 282034, 282035, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282037, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282092, 282708, 282712, 282730, 282730, 282730, 282731, 282731, 282731, 282731, 282732, 282732, 282733, 282733, 282733, 282733, 282734, 282734, 282734, 282734, 282734, 282734, 282734, 282734, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282735, 282800, 282800, 282800, 282800, 282800, 282800, 282800, 282800, 282800, 282800, 282807, 282807, 282807, 282807, 282807, 282807, 282807, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282814, 282842, 282842, 282842};