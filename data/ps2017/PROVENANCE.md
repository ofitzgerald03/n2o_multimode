# Petro and Sedwick (2017) digitised curves

Source: Petro, E. M. and Sedwick, R. J. (2017), "Effects of Water-Vapor Propellant on Electrodeless Thruster Performance", J. Propul. Power 33(6), 1410-1417, doi:10.2514/1.B36389.

Digitised with WebPlotDigitizer, calibrated in data units. The raw exports are in `digitised_raw/`; every other file in this folder is written by `scripts/convert_ps2017.py`. Used by `scripts/ps2017_reproduction.py`.

## Axes

| figure | raw subfolder | x | y | axis scales |
|---|---|---|---|---|
| Fig. 1 | `fig 1/` | T_e [eV] | rate constant kappa [m^3/s] | x linear, y log |
| Fig. 2 | `fig2/` | T_e [eV] | n_i-/n_n [-] | x linear, y linear |
| Fig. 3 | `fig3/` | T_e [eV] | positive-ion fraction [-] | x linear, y linear |
| Fig. 5 | `fig 5/` | T_e [eV] | fractional energy loss chi_j [-] | x linear, y linear |
| Fig. 7a | `fig7a/` | T_e [eV] | thrust efficiency eta_T [-] | x linear, y linear |
| Fig. 7b | `fig7b/` | R_m [-] | thrust efficiency eta_T [-] | x log, y linear |

Fig. 6 (E_c) is not digitised: Petro and Sedwick give analytic fits, Eqs. (9) and (10).

## Raw-to-output map

| raw file | output file | curve | raw rows | output rows | digitising error | flags | sha256 (raw) |
|---|---|---|---|---|---|---|---|
| `fig 1/H2Oplus.csv` | `PS2017_fig1_k_H2Oplus.csv` | H2O+ (ionisation) | 73 | 73 | 0.0126 decades (about 2.9 % relative) | 1 | `c6f0329367e44873634b3d892f601ff0129cb4fac50e4b3773dffa08e491ceda` |
| `fig 1/H2plus.csv` | `PS2017_fig1_k_H2plus.csv` | H2+ (dissociative ionisation) | 57 | 57 | 0.00486 decades (about 1.1 % relative) | 0 | `263f874f3f92de9f8fbd2bbda3e9c28da3bb440672ea27977f1501267543545b` |
| `fig 1/Hplus.csv` | `PS2017_fig1_k_Hplus.csv` | H+ (dissociative ionisation) | 68 | 68 | 0.0132 decades (about 3.1 % relative) | 1 | `4094d071b2b781ed95cc90f99da94e5a2a60142feb67c79aab8f74e3093eaa1b` |
| `fig 1/OHplus.csv` | `PS2017_fig1_k_OHplus.csv` | OH+ (dissociative ionisation) | 59 | 59 | 0.0117 decades (about 2.7 % relative) | 1 | `4e01207babaeaab2c672f462abae18626bacd5cd65cc4a520d071d9efdcca9a6` |
| `fig 1/Oplus.csv` | `PS2017_fig1_k_Oplus.csv` | O+ (dissociative ionisation) | 63 | 63 | 0.0143 decades (about 3.4 % relative) | 0 | `630f0e70026ba25cd190d22a11652f4f15a656a0e098b41b1cdd62c4324f88f3` |
| `fig 1/Oplusplus.csv` | `PS2017_fig1_k_Oplusplus.csv` | O++ (dissociative ionisation) | 21 | 21 | 0.00416 decades (about 0.96 % relative) | 0 | `204155090c41b3f3aced7596ac438a80393b69c47c6655debfa830f61e67d8cc` |
| `fig 1/Hminus.csv` | `PS2017_fig1_k_Hminus.csv` | H- (dissociative attachment) | 61 | 61 | 0.00769 decades (about 1.8 % relative) | 1 | `b69e358a24d41b511d9ec356c7bc37caa13bee43afd6e4dd207c7ceb6f8502e7` |
| `fig 1/OHminus.csv` | `PS2017_fig1_k_OHminus.csv` | OH- (dissociative attachment) | 63 | 63 | 0.00596 decades (about 1.4 % relative) | 1 | `bfd4c7189cdece82b0b0d9147723829ca4c28d25dace99df39acd14731375866` |
| `fig 1/Ominus.csv` | `PS2017_fig1_k_Ominus.csv` | O- (dissociative attachment) | 59 | 59 | 0.00936 decades (about 2.2 % relative) | 1 | `2447a8207994ee1e66489ecab958fa184265f6b3166ce45493bedfc904adde42` |
| `fig2/Petro2017_fig2_line.csv` | `PS2017_fig2_negion_total.csv` | total negative-ion-to-neutral ratio (the single published curve) | 68 | 68 | 7.74e-07 (absolute, -) | 2 | `d62aeef871a16ece87717b3a7908ff85ee5c170da726a3885aed2b4792685bf6` |
| `fig3/H2Oplus.csv` | `PS2017_fig3_frac_H2Oplus.csv` | H2O+ (main axes) | 54 | 54 | 0.000458 (absolute, -) | 0 | `1d3abfbc57a79d3943d185ead4ae41db7e8e79f0d71b47c0aa85beb8a80b10db` |
| `fig3/Hplus.csv` | `PS2017_fig3_frac_Hplus.csv` | H+ (inset) | 44 | 44 | 1.5e-05 (absolute, -) | 0 | `d449564466b08c195f56223c5300e13a060b8e1f29eb7d8439ff2513e99292fb` |
| `fig3/OHplus.csv` | `PS2017_fig3_frac_OHplus.csv` | OH+ (inset) | 41 | 41 | 1.68e-05 (absolute, -) | 0 | `d4f097f1c0beb67a86526fce590bef1745f54e8c7f5194cf832227f44a7c7e91` |
| `fig3/Oplus.csv` | `PS2017_fig3_frac_Oplus.csv` | O+ (inset) | 47 | 47 | 1.4e-05 (absolute, -) | 0 | `da8ea0b00cea5adf37bb6cb3a7880d6f6fa91ef2a646c652c7c6a4cb1ab8f98c` |
| `fig 5/ionization.csv` | `PS2017_fig5_chi_ionisation.csv` | ionisation | 69 | 69 | 0.000404 (absolute, -) | 2 | `c5207e4ab7b4452a9914b3e6fb15642421139363e22ff6d81acef38c09fd2dad` |
| `fig 5/electronic_excitation.csv` | `PS2017_fig5_chi_electronic_excitation.csv` | electronic excitation (incl. radiating and nonradiating dissociation) | 76 | 76 | 0.000654 (absolute, -) | 1 | `608d41e5391e0e47fe025c97a427c6f6da34df9a68481e48ff3583269f32d67a` |
| `fig 5/vibration.csv` | `PS2017_fig5_chi_vibration.csv` | vibration | 70 | 70 | 0.00126 (absolute, -) | 1 | `7818b6463f7079ed4ac66abeac64032963b810842c183e75927db1a4ffe5943e` |
| `fig 5/rotation.csv` | `PS2017_fig5_chi_rotation.csv` | rotation | 78 | 77 | 0.000749 (absolute, -) | 2 | `23430c4fc00e616747f3009ec3c6be0bf0776b5a201688cf81b05b7ef7fc8e1b` |
| `fig 5/elastic.csv` | `PS2017_fig5_chi_elastic.csv` | elastic | 69 | 69 | 0.000633 (absolute, -) | 0 | `41c93e4738103c4e57f2f3f320e772ad0455bc370a94431869c9d8045e2a3b15` |
| `fig 5/attachment.csv` | `PS2017_fig5_chi_attachment.csv` | attachment | 68 | 68 | 0.000501 (absolute, -) | 0 | `4598b989db2d728abc4ac6b7056b1da8eff9f478899c1d7035f908d6da4cc6c9` |
| `fig7a/Ar Rm_1.csv` | `PS2017_fig7a_etaT_Ar_Rm1.csv` | Ar, R_m = 1 | 75 | 75 | 0.000477 (absolute, -) | 0 | `0962d8aaae4ce844456d82f11304a762db934d9a2579f8436e8ca28a0d01396c` |
| `fig7a/Ar Rm_10.csv` | `PS2017_fig7a_etaT_Ar_Rm10.csv` | Ar, R_m = 10 | 75 | 75 | 0.000787 (absolute, -) | 0 | `2f2417c61f1fd48c3395c1b94d1a49ece6a0196831e77b0c918c8a98662cdacc` |
| `fig7a/Ar Rm_100.csv` | `PS2017_fig7a_etaT_Ar_Rm100.csv` | Ar, R_m = 100 | 74 | 74 | 0.000541 (absolute, -) | 0 | `ee465b985f1ec885097605f5f178118d9b27eaeece180f91588e2a5c486cb603` |
| `fig7a/H2O Rm_1.csv` | `PS2017_fig7a_etaT_H2O_Rm1.csv` | H2O, R_m = 1 | 68 | 68 | 0.00119 (absolute, -) | 1 | `543c003173cae7a306c8be157a4e5a4ad4427ac131cf2511c9a945289294d321` |
| `fig7a/H2O Rm_10.csv` | `PS2017_fig7a_etaT_H2O_Rm10.csv` | H2O, R_m = 10 | 43 | 43 | 0.00111 (absolute, -) | 0 | `a092d22d02dab321e647b12a1819ddb67782c821b3f75996d7ab9ca8a6128dca` |
| `fig7a/H2O Rm_100.csv` | `PS2017_fig7a_etaT_H2O_Rm100.csv` | H2O, R_m = 100 | 43 | 43 | 0.00126 (absolute, -) | 0 | `fc5983ec320067f31dcfcb52936c4598bb27cfa74e5a5499c45f9ee5a4a0fc1a` |
| `fig7b/Ar Te_2.csv` | `PS2017_fig7b_etaT_Ar_Te2.csv` | Ar, T_e = 2 eV | 73 | 73 | 0.000233 (absolute, -) | 0 | `3bc25bd616cc6699bdd923d86149593b7a984a204c4843b69559646cfb6bd4e3` |
| `fig7b/Ar Te_5.csv` | `PS2017_fig7b_etaT_Ar_Te5.csv` | Ar, T_e = 5 eV | 73 | 73 | 0.000385 (absolute, -) | 0 | `1e0712f3df30e589e96515ba3d352583e89e377283803090a75a3a5e8a7b34dd` |
| `fig7b/Ar Te_10.csv` | `PS2017_fig7b_etaT_Ar_Te10.csv` | Ar, T_e = 10 eV | 73 | 73 | 0.00062 (absolute, -) | 0 | `db17099c81327445419ba305254a2e6a1806cb36a27b533c48ad13c03e189710` |
| `fig7b/H2O Te_2.csv` | `PS2017_fig7b_etaT_H2O_Te2.csv` | H2O, T_e = 2 eV | 40 | 40 | 0.000533 (absolute, -) | 0 | `151ca85b10bf2e30f929030cb39e94b5c76c9998687e3bda91fd1b06893bcfce` |
| `fig7b/H2O Te_5.csv` | `PS2017_fig7b_etaT_H2O_Te5.csv` | H2O, T_e = 5 eV | 40 | 40 | 0.000737 (absolute, -) | 0 | `1cfc599293d538be54ce0d7b055159799d42d8a1d074264cc8490f729749a870` |
| `fig7b/H2O Te_10.csv` | `PS2017_fig7b_etaT_H2O_Te10.csv` | H2O, T_e = 10 eV | 40 | 40 | 0.000563 (absolute, -) | 0 | `a63379b2e0ae26eb61fe93098c65f36d3c386d4e369db84e23b33047fba44107` |

## Removed points

- `PS2017_fig5_chi_rotation`: removed (4.286632014291518, 0.05737715325280757): stray click (raw row 35; 0.013 below both neighbours)

## Flags (kept as digitised unless stated)

- `PS2017_fig1_k_H2Oplus`:
    - residual outlier at x = 1.6172 (y = 2.72027e-18): +5.9 x RMS
- `PS2017_fig1_k_Hplus`:
    - residual outlier at x = 1.62626 (y = 4.02345e-21): +5.5 x RMS
- `PS2017_fig1_k_OHplus`:
    - residual outlier at x = 2.32365 (y = 1.08517e-20): +4.2 x RMS
- `PS2017_fig1_k_Hminus`:
    - residual outlier at x = 1.30726 (y = 1.14418e-17): -4.0 x RMS
- `PS2017_fig1_k_OHminus`:
    - residual outlier at x = 1.67689 (y = 8.39301e-19): +4.9 x RMS
- `PS2017_fig1_k_Ominus`:
    - residual outlier at x = 1.61073 (y = 2.01885e-18): +4.5 x RMS
- `PS2017_fig2_negion_total`:
    - residual outlier at x = 1.09734 (y = 0.000566553): +4.3 x RMS
    - residual outlier at x = 1.37229 (y = 0.00066204): +4.3 x RMS
- `PS2017_fig5_chi_ionisation`:
    - residual outlier at x = 1.66003 (y = 0.00636649): -4.2 x RMS
    - residual outlier at x = 3.19037 (y = 0.193615): +4.1 x RMS
- `PS2017_fig5_chi_electronic_excitation`:
    - residual outlier at x = 2.31262 (y = 0.297108): +5.9 x RMS
- `PS2017_fig5_chi_vibration`:
    - residual outlier at x = 2.88222 (y = 0.206886): -5.1 x RMS
- `PS2017_fig5_chi_rotation`:
    - residual outlier at x = 2.57503 (y = 0.266112): +4.8 x RMS
    - residual outlier at x = 2.73628 (y = 0.227328): -5.0 x RMS
- `PS2017_fig7a_etaT_H2O_Rm1`:
    - negative y at x = 1.1653744311860157 (y = -0.0005919921145212115): clamped to 0

## Fridman detachment curves

Electron-impact detachment cross sections of H- and O-, Fridman (2021), Fig. 2.9 (cited by Petro and Sedwick for detachment), digitised; in `data/cross_sections/n2o/digitised_raw/Fridman2021_fig2_9/`, two columns without a header, E [eV] and sigma [cm^2]. Used by `scripts/ps2017_reproduction.py`.
