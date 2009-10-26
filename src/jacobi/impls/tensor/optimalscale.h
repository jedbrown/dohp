/**
* @file   optimalscale.h
* @author Jed Brown <jed@59A2.org>
* @date   Sun May 10 14:18:01 2009
* 
* @brief  Defines optimal scaling for q1 matrices assembled on Legendre-Gauss_Lobatto points.
*
* These arrays were generated by optimizing over all scalings of the form
*
* As = diag(s) * A * diag(s)
*
* where A is either the Q1 mass or Laplacian matrix, to minimize the "effective condition number" of inv(As)*A.  The
* "effective condition number" is computed from the SVD by dropping singular values that are 0 (and would be supplied by
* boundary values).
*
**/

#if !defined(_OPTIMALSCALE_H)
#define _OPTIMALSCALE_H

#include <dohptype.h>

static const dReal optimal_ones[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

static const dReal optimal_mscale_2[] = {1,1},
  optimal_mscale_3[] = {0.8633789062499998, 1.2732421875000004, 0.8633789062499998},
  optimal_mscale_4[] = {0.87558593749999991, 1.1244140625000001, 1.1244140625000001, 0.87558593749999991},
  optimal_mscale_5[] = {0.84529011510312568, 1.1746890407055619, 0.96004168838262516, 1.1746890407055619, 0.84529011510312568},
  optimal_mscale_6[] = {0.83414557335072448, 1.1347409904428778, 1.0311134362063976, 1.0311134362063976, 1.1347409904428778, 0.83414557335072448},
  optimal_mscale_7[] = {0.78555955394946975, 1.1276363257065061, 1.0360325284064045, 1.1015431838752399, 1.0360325284064045, 1.1276363257065061, 0.78555955394946975},
  optimal_mscale_8[] = {0.78005408267437315, 1.115947360985226, 1.0720759672236908, 1.0319225891167099, 1.0319225891167099, 1.0720759672236908, 1.115947360985226, 0.78005408267437315},
  optimal_mscale_9[] = {0.76452030781575819, 1.0943844400157769, 1.065753105798533, 1.1126528480337201, 0.92537859667242406, 1.1126528480337201, 1.065753105798533, 1.0943844400157769, 0.76452030781575819},
  optimal_mscale_10[] = {0.76899262661096279, 1.0844808331846316, 1.0650732246613432, 1.1092630792818661, 0.97219023626119672, 0.97219023626119672, 1.1092630792818661, 1.0650732246613432, 1.0844808331846316, 0.76899262661096279},
  optimal_mscale_11[] = {0.7359699636059529, 1.0641021099995149, 1.0486739809231413, 1.0951863533499089, 1.0124160181349324, 1.0873031479730995, 1.0124160181349324, 1.0951863533499089, 1.0486739809231413, 1.0641021099995149, 0.7359699636059529},
  optimal_mscale_12[] = {0.75147025229575348, 1.0656982122596839, 1.0354804945478802, 1.0969247634799855, 1.0244541215961656, 1.0259721558205308, 1.0259721558205308, 1.0244541215961656, 1.0969247634799855, 1.0354804945478802, 1.0656982122596839, 0.75147025229575348},
  optimal_mscale_13[] = {0.73513052753657493, 1.0549675638391163, 1.0187361597988762, 1.0881864667657033, 1.0221852293656799, 1.0598735701897892, 1.0418409650085196, 1.0598735701897892, 1.0221852293656799, 1.0881864667657033, 1.0187361597988762, 1.0549675638391163, 0.73513052753657493},
  optimal_mscale_14[] = {0.74324111138944104, 1.0579035612110257, 0.99754899979983325, 1.0843431193507969, 1.0105256581725728, 1.0340452975350765, 1.0723922525412544, 1.0723922525412544, 1.0340452975350765, 1.0105256581725728, 1.0843431193507969, 0.99754899979983325, 1.0579035612110257, 0.74324111138944104};

static const dReal optimal_lscale_2[] = {1,1},
  optimal_lscale_3[] = {0.93129882812499987, 1.1374023437500003, 0.93129882812499987},
  optimal_lscale_4[] = {0.95996093749999978, 1.0400390625000002, 1.0400390625000002, 0.95996093749999978},
  optimal_lscale_5[] = {0.95323834643054373, 1.0013589126035178, 1.0908054819318769, 1.0013589126035178, 0.95323834643054373},
  optimal_lscale_6[] = {0.97233762780379052, 0.99478205400373554, 1.0328803181924742, 1.0328803181924742, 0.99478205400373554, 0.97233762780379052},
  optimal_lscale_7[] = {0.97491271419538439, 0.9875817846736441, 1.0112409122406048, 1.0525291777807337, 1.0112409122406048, 0.9875817846736441, 0.97491271419538439},
  optimal_lscale_8[] = {0.98242904574845613, 0.98945640534627399, 1.0008163977013877, 1.0272981512038823, 1.0272981512038823, 1.0008163977013877, 0.98945640534627399, 0.98242904574845613},
  optimal_lscale_9[] = {0.98371234834817733, 0.98934385044409012, 0.99851920367213542, 1.0116345126548261, 1.0335801697615423, 1.0116345126548261, 0.99851920367213542, 0.98934385044409012, 0.98371234834817733},
  optimal_lscale_10[] = {0.98713874330193863, 0.99110664000610571, 0.99544593877259457, 1.0054188667757713, 1.02088981114359, 1.02088981114359, 1.0054188667757713, 0.99544593877259457, 0.99110664000610571, 0.98713874330193863},
  optimal_lscale_11[] = {0.98802560946845452, 0.99088308151773186, 0.99485717770327386, 1.0018012691703861, 1.0132746372828709, 1.022316449714566, 1.0132746372828709, 1.0018012691703861, 0.99485717770327386, 0.99088308151773186, 0.98802560946845452},
  optimal_lscale_12[] = {0.98933384410274483, 0.99148492407376065, 0.99456924585078266, 0.99974133924938646, 1.0076339300810422, 1.0172367166422831, 1.0172367166422831, 1.0076339300810422, 0.99974133924938646, 0.99456924585078266, 0.99148492407376065, 0.98933384410274483},
  optimal_lscale_13[] = {0.99268113885731557, 0.99445091830329768, 0.99613360268945739, 0.99536834490944437, 1.0054811657260543, 1.0081328020609224, 1.0155040549070158, 1.0081328020609224, 1.0054811657260543, 0.99536834490944437, 0.99613360268945739, 0.99445091830329768, 0.99268113885731557},
  optimal_lscale_14[] = {0.99091904296866629, 0.99240628966641609, 0.99441004182561787, 0.99649973126893876, 1.0006824544025643, 1.0087390373686946, 1.016343402499102, 1.016343402499102, 1.0087390373686946, 1.0006824544025643, 0.99649973126893876, 0.99441004182561787, 0.99240628966641609, 0.99091904296866629};

#endif
