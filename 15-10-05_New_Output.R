### Metabolites with more than 50% missing values for each dataset ###

> metab.nas(SynExp1)
[1] "Arginine"
> metab.nas(SynExp2)
[1] "Alanine"               "Arginine"              "Glucose-6-phosphate"  
[4] "Glutaric acid, 2-oxo-" "Isoleucine"            "Ornithine"            
[7] "Pyruvic acid"          "Serine"                "Succinic acid"        
[10] "Threonine"             "Tyrosine"              "Valine"               
> metab.nas(SynExp3)
[1] "Glyceric acid" "Pyruvic acid"  "Succinic acid" "Valine"       
> metab.nas(SynExp4)
[1] "Lysine"
> metab.nas(AraExp1day)
[1] "Glyceric acid-2-phosphate" "Glyceric acid-3-phosphate" "Glycolic acid-2-phosphate"
[4] "Phosphoenolpyruvic acid"  
> metab.nas(AraExp1night)
[1] "Glyceric acid-2-phosphate" "Glyceric acid-3-phosphate" "Glycolic acid-2-phosphate"
[4] "Phosphoenolpyruvic acid"  
> metab.nas(AraExp2VC)
[1] "Arginine"                  "Benzoic acid"              "Ethanolamine"             
[4] "Glyceric acid-2-phosphate" "Glyceric acid-3-phosphate" "Glycerol"                 
[7] "Glycolic acid-2-phosphate" "Lysine"                    "Phosphoenolpyruvic acid"  
> metab.nas(AraExp2HC)
[1] "Arginine"                  "Benzoic acid"              "Ethanolamine"             
[4] "Glyceric acid-2-phosphate" "Glyceric acid-3-phosphate" "Glycerol"                 
[7] "Glycolic acid-2-phosphate" "Lysine" 

### Percentage of missing values after removing metabolites with more than 50% missing values ###

> percent.nas(SynExp1.1)
[1] 0
> percent.nas(SynExp3.1)
[1] 0.4705882
> percent.nas(SynExp4.1)
[1] 0
> percent.nas(AraExp1day.1)
[1] 0.9803922
> percent.nas(AraExp1night.1)
[1] 3.431373
> percent.nas(AraExp2VC.1)
[1] 1.641587
> percent.nas(AraExp2HC.1)
[1] 4.061625

### Average silhouette coefficient to choose cluster number for PAM biclustering ###

> num.clus(Scaled)
2  avg:  0.07 clus:  0.04 0.1 
3  avg:  0.08 clus:  0 0.07 0.12 
4  avg:  0.09 clus:  0 0.03 0.2 0 
5  avg:  0.07 clus:  0 -0.03 0 0.28 0 
6  avg:  0.07 clus:  0 -0.02 0 0.27 0 0 
7  avg:  0.07 clus:  0 -0.01 0 0 0.26 0 0 
8  avg:  0.07 clus:  0 0 0 0 0.24 0 0 0 
9  avg:  0.08 clus:  0 0.05 0 0 0 0.22 0 0 0 
10  avg:  0.09 clus:  0 0.12 0 0 0 0.2 0 0 0 0 
> num.clus(t(Scaled))
2  avg:  0.15 clus:  0.14 0.16 
3  avg:  0.14 clus:  0.11 0.16 0.17 
4  avg:  0.15 clus:  0.1 0.14 0.16 0.18 
5  avg:  0.14 clus:  0.09 0.12 0.16 0.14 0.2 
6  avg:  0.15 clus:  0.16 0.13 0.18 0.09 0.08 0.21 
7  avg:  0.15 clus:  0.2 0.12 0.14 0.08 0.07 0.23 0.21 
8  avg:  0.15 clus:  0.2 0.12 0.14 0.08 0.08 0.1 0.2 0.25 
9  avg:  0.14 clus:  0.05 0.2 0.13 0.08 0.17 0.08 0.08 0.2 0.25 
10  avg:  0.14 clus:  0.04 0.32 0.11 0.13 0.16 0.07 0.07 0.19 0.28 0.2 
> num.clus(Transformed)
2  avg:  0.37 clus:  0.39 0 
3  avg:  0.05 clus:  0.03 0.1 0 
4  avg:  0.08 clus:  0 0.05 0.11 0 
5  avg:  0.07 clus:  0 0.12 0.06 0 0 
6  avg:  0.07 clus:  0 0.07 0.12 0 0 0 
7  avg:  0.06 clus:  0 0.1 0.07 0 0 0 0 
8  avg:  0.07 clus:  0 0.08 0.13 0 0 0 0 0 
9  avg:  0.07 clus:  0 0.05 0 0 0 0.2 0 0 0 
10  avg:  0.06 clus:  0 0.05 0 0 0 0 0.17 0 0 0 
> num.clus(t(Transformed))
2  avg:  0.4 clus:  0.4 0.42 
3  avg:  0.19 clus:  0.21 0.12 0.45 
4  avg:  0.16 clus:  0.14 0.08 0.17 0.41 
5  avg:  0.16 clus:  0.13 0.13 0.15 0.12 0.4 
6  avg:  0.17 clus:  0.39 0.12 0.14 0.12 0.16 0.4 
7  avg:  0.16 clus:  0.37 0.11 0.12 0.11 0.11 0.4 0.25 
8  avg:  0.16 clus:  0.37 0.06 0.12 0.2 0.1 0.38 0.11 0.29 
9  avg:  0.16 clus:  0.36 0.07 0.11 0.17 0.09 0.18 0.08 0.29 0.41 
10  avg:  0.16 clus:  0.36 0.06 0.09 0.17 0.09 0.18 0.07 0.28 0.34 0.5

### PAM clusters of samples ###

> S.clu1 <- S.PAM[S.PAM[,ncol(S.PAM)]==1,]
> rownames(S.clu1) 
[1] "X10076if_8_WT.HC"           "X10076if_17_WT.HC"         
[3] "X10076if_25_WT.HC"          "X10076if_34_WT.HC"         
[5] "X10076if_20_WT.3hLC"        "X10076if_28_WT.3hLC"       
[7] "X10076if_37_WT.3hLC"        "X13122mo_6_WT_A_P1_0h_LC"  
[9] "X13122mo_66_WT_A_P2_0h_LC"  "X13122mo_30_WT_A_P3_0h_LC" 
[11] "X13122mo_42_WT_B_P1_0h_LC"  "X13122mo_54_WT_B_P2_0h_LC" 
[13] "X13127mo_10_WT_B_P3_0h_LC"  "X13127mo_21_WT_C_P1_0h_LC" 
[15] "X13127mo_32_WT_C_P2_0h_LC"  "X13127mo_44_WT_C_P3_0h_LC" 
[17] "X13122mo_58_WT_B_P2_3h_LC"  "X13122mo_62_WT_A_P1_24h_LC"
[19] "X13122mo_38_WT_A_P3_24h_LC" "X10076if_8_WT.HC"          
[21] "X10076if_17_WT.HC"          "X10076if_34_WT.HC"         
[23] "X10076if_20_WT.3hLC"        "X10076if_28_WT.3hLC"       
[25] "X10076if_37_WT.3hLC"        "X10152oa_3_col.0.HC.14.1"  
[27] "X10152oa_24_col.0.HC.14.2"  "X10152oa_45_col.0.HC.14.3" 
[29] "X10152oa_25_col.0.LC1.14.2" "X10152oa_46_col.0.LC1.14.3"
[31] "X10153oa_25_col.0.LC5.14.2" "X10152oa_26_col.0.HC.20.2" 
[33] "X10152oa_6_col.0.LC1.20.1"  "X10152oa_27_col.0.LC1.20.2"
[35] "X10152oa_48_col.0.LC1.20.3" "X10153oa_5_col.0.LC3.20.1" 
[37] "X10153oa_26_col.0.LC3.20.2" "X10153oa_47_col.0.LC3.20.3"
[39] "X10153oa_6_col.0.LC5.20.1"  "X10153oa_48_col.0.LC5.20.3"
[41] "X10201hl_13_NC6.1"          "X10201hl_14_VC6.1"         
[43] "X10201hl_15_NC7.1"          "X10201hl_16_VC7.1"         
[45] "X10201hl_19_VC1.2"          "X10201hl_21_VC2.2"         
[47] "X10201hl_29_VC6.2"          "X10201hl_30_NC7.2"         
[49] "X10201hl_31_VC7.2"          "X10201hl_33_NC1.3"         
[51] "X10201hl_34_VC1.3"          "X10201hl_36_VC2.3"         
[53] "X10201hl_38_VC3.3"          "X10201hl_4_VC1.1"          
[55] "X10201hl_44_VC6.3"          "X10201hl_45_NC7.3"         
[57] "X10201hl_46_VC7.3"          "X10201hl_48_VC3.3.2"       
[59] "X10201hl_6_VC2.1"           "X10201hl_8_VC3.1"          
[61] "X10202hl_14_VC6.4"          "X10202hl_15_NC7.4"         
[63] "X10202hl_16_VC7.4"          "X10202hl_18_NC1.5"         
[65] "X10202hl_19_VC1.5"          "X10202hl_21_VC2.5"         
[67] "X10202hl_29_VC6.5"          "X10202hl_30_NC7.5"         
[69] "X10202hl_31_VC7.5"          "X10202hl_33_NC1.6"         
[71] "X10202hl_34_VC1.6"          "X10202hl_36_VC2.6"         
[73] "X10202hl_38_VC3.6"          "X10202hl_4_VC1.4"          
[75] "X10202hl_43_NC6.6"          "X10202hl_44_VC6.6"         
[77] "X10202hl_45_NC7.6"          "X10202hl_46_VC7.6"         
[79] "X10202hl_48_VC3.6.2"        "X10202hl_6_VC2.4"          
[81] "X10202hl_8_VC3.4"           "X10193hl_15_HC6.1"         
[83] "X10193hl_16_NC7.1"          "X10193hl_17_HC7.1"         
[85] "X10193hl_19_HC1.2"          "X10193hl_21_HC2.2"         
[87] "X10193hl_29_NC6.2"          "X10193hl_30_HC6.2"         
[89] "X10193hl_31_NC7.2"          "X10193hl_32_HC7.2"         
[91] "X10193hl_33_NC1.3"          "X10193hl_34_HC1.3"         
[93] "X10193hl_36_HC2.3"          "X10193hl_4_HC1.1"          
[95] "X10193hl_41_HC4.3"          "X10193hl_44_NC6.3"         
[97] "X10193hl_45_HC6.3"          "X10193hl_46_HC7.3"         
[99] "X10193hl_47_NC7.3"          "X10194hl_15_HC6.4"         
[101] "X10194hl_16_NC7.4"          "X10194hl_17_HC7.4"         
[103] "X10194hl_19_HC1.5"          "X10194hl_21_HC2.5"         
[105] "X10194hl_29_NC6.5"          "X10194hl_30_HC6.5"         
[107] "X10194hl_31_NC7.5"          "X10194hl_32_HC7.5"         
[109] "X10194hl_33_NC1.6"          "X10194hl_34_HC1.6"         
[111] "X10194hl_4_HC1.4"           "X10194hl_44_NC6.6"         
[113] "X10194hl_45_HC6.6"          "X10194hl_46_HC7.6"         
[115] "X10194hl_47_NC7.6"         
> 
  > Tr.clu1 <- Tr.PAM[Tr.PAM[,ncol(Tr.PAM)]==1,]
> rownames(Tr.clu1)
[1] "X10076if_8_WT.HC"           "X10076if_17_WT.HC"         
[3] "X10076if_25_WT.HC"          "X10076if_34_WT.HC"         
[5] "X10076if_11_WT.3hLC"        "X10076if_20_WT.3hLC"       
[7] "X10076if_28_WT.3hLC"        "X10076if_37_WT.3hLC"       
[9] "X10076if_14_WT.24hLC"       "X10076if_23_WT.24hLC"      
[11] "X10076if_31_WT.24hLC"       "X10076if_40_WT.24hLC"      
[13] "X13122mo_6_WT_A_P1_0h_LC"   "X13122mo_66_WT_A_P2_0h_LC" 
[15] "X13122mo_30_WT_A_P3_0h_LC"  "X13122mo_42_WT_B_P1_0h_LC" 
[17] "X13122mo_54_WT_B_P2_0h_LC"  "X13127mo_10_WT_B_P3_0h_LC" 
[19] "X13127mo_21_WT_C_P1_0h_LC"  "X13127mo_32_WT_C_P2_0h_LC" 
[21] "X13127mo_44_WT_C_P3_0h_LC"  "X13122mo_10_WT_A_P1_3h_LC" 
[23] "X13122mo_22_WT_A_P2_3h_LC"  "X13122mo_34_WT_A_P3_3h_LC" 
[25] "X13122mo_46_WT_B_P1_3h_LC"  "X13122mo_58_WT_B_P2_3h_LC" 
[27] "X13127mo_36_WT_C_P2_3h_LC"  "X13127mo_48_WT_C_P3_3h_LC" 
[29] "X13122mo_26_WT_A_P2_24h_LC" "X13122mo_38_WT_A_P3_24h_LC"
[31] "X13122mo_50_WT_B_P1_24h_LC" "X13127mo_6_WT_B_P2_24h_LC" 
[33] "X13127mo_17_WT_B_P3_24h_LC" "X13127mo_28_WT_C_P1_24h_LC"
[35] "X13127mo_40_WT_C_P2_24h_LC" "X13127mo_52_WT_C_P3_24h_LC"
[37] "X10076if_8_WT.HC"           "X10076if_17_WT.HC"         
[39] "X10076if_25_WT.HC"          "X10076if_34_WT.HC"         
[41] "X10076if_11_WT.3hLC"        "X10076if_20_WT.3hLC"       
[43] "X10076if_28_WT.3hLC"        "X10076if_37_WT.3hLC"       
[45] "X10076if_14_WT.24hLC"       "X10076if_23_WT.24hLC"      
[47] "X10152oa_3_col.0.HC.14.1"   "X10152oa_24_col.0.HC.14.2" 
[49] "X10152oa_45_col.0.HC.14.3"  "X10152oa_4_col.0.LC1.14.1" 
[51] "X10152oa_25_col.0.LC1.14.2" "X10152oa_46_col.0.LC1.14.3"
[53] "X10153oa_3_col.0.LC3.14.1"  "X10153oa_24_col.0.LC3.14.2"
[55] "X10153oa_45_col.0.LC3.14.3" "X10153oa_4_col.0.LC5.14.1" 
[57] "X10153oa_25_col.0.LC5.14.2" "X10153oa_46_col.0.LC5.14.3"
[59] "X10152oa_5_col.0.HC.20.1"   "X10152oa_26_col.0.HC.20.2" 
[61] "X10152oa_47_col.0.HC.20.3"  "X10152oa_6_col.0.LC1.20.1" 
[63] "X10152oa_27_col.0.LC1.20.2" "X10152oa_48_col.0.LC1.20.3"
[65] "X10153oa_5_col.0.LC3.20.1"  "X10153oa_26_col.0.LC3.20.2"
[67] "X10153oa_47_col.0.LC3.20.3" "X10153oa_6_col.0.LC5.20.1" 
[69] "X10153oa_27_col.0.LC5.20.2" "X10153oa_48_col.0.LC5.20.3"
[71] "X10201hl_10_VC4.1"          "X10201hl_11_NC5.1"         
[73] "X10201hl_12_VC5.1"          "X10201hl_13_NC6.1"         
[75] "X10201hl_14_VC6.1"          "X10201hl_18_NC1.2"         
[77] "X10201hl_19_VC1.2"          "X10201hl_20_NC2.2"         
[79] "X10201hl_21_VC2.2"          "X10201hl_22_NC3.2"         
[81] "X10201hl_23_VC3.2"          "X10201hl_24_NC4.2"         
[83] "X10201hl_25_VC4.2"          "X10201hl_26_NC5.2"         
[85] "X10201hl_27_VC5.2"          "X10201hl_28_NC6.2"         
[87] "X10201hl_29_VC6.2"          "X10201hl_3_NC1.1"          
[89] "X10201hl_33_NC1.3"          "X10201hl_34_VC1.3"         
[91] "X10201hl_35_NC2.3"          "X10201hl_36_VC2.3"         
[93] "X10201hl_37_NC3.3"          "X10201hl_38_VC3.3"         
[95] "X10201hl_39_NC4.3"          "X10201hl_4_VC1.1"          
[97] "X10201hl_40_VC4.3"          "X10201hl_41_NC5.3"         
[99] "X10201hl_42_VC5.3"          "X10201hl_43_NC6.3"         
[101] "X10201hl_44_VC6.3"          "X10201hl_48_VC3.3.2"       
[103] "X10201hl_5_NC2.1"           "X10201hl_6_VC2.1"          
[105] "X10201hl_7_NC3.1"           "X10201hl_8_VC3.1"          
[107] "X10201hl_9_NC4.1"           "X10202hl_10_VC4.4"         
[109] "X10202hl_11_NC5.4"          "X10202hl_12_VC5.4"         
[111] "X10202hl_13_NC6.4"          "X10202hl_14_VC6.4"         
[113] "X10202hl_18_NC1.5"          "X10202hl_19_VC1.5"         
[115] "X10202hl_20_NC2.5"          "X10202hl_21_VC2.5"         
[117] "X10202hl_22_NC3.5"          "X10202hl_23_VC3.5"         
[119] "X10202hl_24_NC4.5"          "X10202hl_25_VC4.5"         
[121] "X10202hl_26_NC5.5"          "X10202hl_27_VC5.5"         
[123] "X10202hl_28_NC6.5"          "X10202hl_29_VC6.5"         
[125] "X10202hl_3_NC1.4"           "X10202hl_33_NC1.6"         
[127] "X10202hl_34_VC1.6"          "X10202hl_35_NC2.6"         
[129] "X10202hl_36_VC2.6"          "X10202hl_37_NC3.6"         
[131] "X10202hl_38_VC3.6"          "X10202hl_39_NC4.6"         
[133] "X10202hl_4_VC1.4"           "X10202hl_40_VC4.6"         
[135] "X10202hl_41_NC5.6"          "X10202hl_42_VC5.6"         
[137] "X10202hl_43_NC6.6"          "X10202hl_48_VC3.6.2"       
[139] "X10202hl_5_NC2.4"           "X10202hl_6_VC2.4"          
[141] "X10202hl_7_NC3.4"           "X10202hl_8_VC3.4"          
[143] "X10202hl_9_NC4.4"           "X10193hl_10_HC4.1"         
[145] "X10193hl_11_NC5.1"          "X10193hl_12_HC5.1"         
[147] "X10193hl_14_NC6.1"          "X10193hl_18_NC1.2"         
[149] "X10193hl_19_HC1.2"          "X10193hl_20_NC2.2"         
[151] "X10193hl_21_HC2.2"          "X10193hl_22_NC3.2"         
[153] "X10193hl_23_HC3.2"          "X10193hl_24_NC4.2"         
[155] "X10193hl_25_HC4.2"          "X10193hl_27_NC5.2"         
[157] "X10193hl_28_HC5.2"          "X10193hl_29_NC6.2"         
[159] "X10193hl_3_NC1.1"           "X10193hl_30_HC6.2"         
[161] "X10193hl_33_NC1.3"          "X10193hl_34_HC1.3"         
[163] "X10193hl_35_NC2.3"          "X10193hl_36_HC2.3"         
[165] "X10193hl_38_NC3.3"          "X10193hl_39_HC3.3"         
[167] "X10193hl_4_HC1.1"           "X10193hl_40_NC4.3"         
[169] "X10193hl_41_HC4.3"          "X10193hl_42_NC5.3"         
[171] "X10193hl_43_HC5.3"          "X10193hl_44_NC6.3"         
[173] "X10193hl_5_NC2.1"           "X10193hl_6_HC2.1"          
[175] "X10193hl_7_NC3.1"           "X10193hl_8_HC3.1"          
[177] "X10193hl_9_NC4.1"           "X10194hl_10_HC4.4"         
[179] "X10194hl_11_NC5.4"          "X10194hl_12_HC5.4"         
[181] "X10194hl_14_NC6.4"          "X10194hl_15_HC6.4"         
[183] "X10194hl_18_NC1.5"          "X10194hl_19_HC1.5"         
[185] "X10194hl_20_NC2.5"          "X10194hl_21_HC2.5"         
[187] "X10194hl_22_NC3.5"          "X10194hl_23_HC3.5"         
[189] "X10194hl_24_NC4.5"          "X10194hl_25_HC4.5"         
[191] "X10194hl_27_NC5.5"          "X10194hl_28_HC5.5"         
[193] "X10194hl_29_NC6.5"          "X10194hl_3_NC1.4"          
[195] "X10194hl_30_HC6.5"          "X10194hl_33_NC1.6"         
[197] "X10194hl_34_HC1.6"          "X10194hl_35_NC2.6"         
[199] "X10194hl_36_HC2.6"          "X10194hl_38_NC3.6"         
[201] "X10194hl_39_HC3.6"          "X10194hl_4_HC1.4"          
[203] "X10194hl_40_NC4.6"          "X10194hl_41_HC4.6"         
[205] "X10194hl_42_NC5.6"          "X10194hl_43_HC5.6"         
[207] "X10194hl_44_NC6.6"          "X10194hl_45_HC6.6"         
[209] "X10194hl_5_NC2.4"           "X10194hl_6_HC2.4"          
[211] "X10194hl_7_NC3.4"           "X10194hl_8_HC3.4"          
[213] "X10194hl_9_NC4.4"          
> 
  > S.clu2 <- S.PAM[S.PAM[,ncol(S.PAM)]==2,]
> rownames(S.clu2)
[1] "X10076if_11_WT.3hLC"        "X10076if_14_WT.24hLC"      
[3] "X10076if_23_WT.24hLC"       "X10076if_31_WT.24hLC"      
[5] "X10076if_40_WT.24hLC"       "X13122mo_10_WT_A_P1_3h_LC" 
[7] "X13122mo_22_WT_A_P2_3h_LC"  "X13122mo_34_WT_A_P3_3h_LC" 
[9] "X13122mo_46_WT_B_P1_3h_LC"  "X13127mo_36_WT_C_P2_3h_LC" 
[11] "X13127mo_48_WT_C_P3_3h_LC"  "X13122mo_26_WT_A_P2_24h_LC"
[13] "X13122mo_50_WT_B_P1_24h_LC" "X13127mo_6_WT_B_P2_24h_LC" 
[15] "X13127mo_17_WT_B_P3_24h_LC" "X13127mo_28_WT_C_P1_24h_LC"
[17] "X13127mo_40_WT_C_P2_24h_LC" "X13127mo_52_WT_C_P3_24h_LC"
[19] "X10076if_25_WT.HC"          "X10076if_11_WT.3hLC"       
[21] "X10076if_14_WT.24hLC"       "X10076if_23_WT.24hLC"      
[23] "X10152oa_4_col.0.LC1.14.1"  "X10153oa_3_col.0.LC3.14.1" 
[25] "X10153oa_24_col.0.LC3.14.2" "X10153oa_45_col.0.LC3.14.3"
[27] "X10153oa_4_col.0.LC5.14.1"  "X10153oa_46_col.0.LC5.14.3"
[29] "X10152oa_5_col.0.HC.20.1"   "X10152oa_47_col.0.HC.20.3" 
[31] "X10153oa_27_col.0.LC5.20.2" "X10201hl_10_VC4.1"         
[33] "X10201hl_11_NC5.1"          "X10201hl_12_VC5.1"         
[35] "X10201hl_18_NC1.2"          "X10201hl_20_NC2.2"         
[37] "X10201hl_22_NC3.2"          "X10201hl_23_VC3.2"         
[39] "X10201hl_24_NC4.2"          "X10201hl_25_VC4.2"         
[41] "X10201hl_26_NC5.2"          "X10201hl_27_VC5.2"         
[43] "X10201hl_28_NC6.2"          "X10201hl_3_NC1.1"          
[45] "X10201hl_35_NC2.3"          "X10201hl_37_NC3.3"         
[47] "X10201hl_39_NC4.3"          "X10201hl_40_VC4.3"         
[49] "X10201hl_41_NC5.3"          "X10201hl_42_VC5.3"         
[51] "X10201hl_43_NC6.3"          "X10201hl_5_NC2.1"          
[53] "X10201hl_7_NC3.1"           "X10201hl_9_NC4.1"          
[55] "X10202hl_10_VC4.4"          "X10202hl_11_NC5.4"         
[57] "X10202hl_12_VC5.4"          "X10202hl_13_NC6.4"         
[59] "X10202hl_20_NC2.5"          "X10202hl_22_NC3.5"         
[61] "X10202hl_23_VC3.5"          "X10202hl_24_NC4.5"         
[63] "X10202hl_25_VC4.5"          "X10202hl_26_NC5.5"         
[65] "X10202hl_27_VC5.5"          "X10202hl_28_NC6.5"         
[67] "X10202hl_3_NC1.4"           "X10202hl_35_NC2.6"         
[69] "X10202hl_37_NC3.6"          "X10202hl_39_NC4.6"         
[71] "X10202hl_40_VC4.6"          "X10202hl_41_NC5.6"         
[73] "X10202hl_42_VC5.6"          "X10202hl_5_NC2.4"          
[75] "X10202hl_7_NC3.4"           "X10202hl_9_NC4.4"          
[77] "X10193hl_10_HC4.1"          "X10193hl_11_NC5.1"         
[79] "X10193hl_12_HC5.1"          "X10193hl_14_NC6.1"         
[81] "X10193hl_18_NC1.2"          "X10193hl_20_NC2.2"         
[83] "X10193hl_22_NC3.2"          "X10193hl_23_HC3.2"         
[85] "X10193hl_24_NC4.2"          "X10193hl_25_HC4.2"         
[87] "X10193hl_27_NC5.2"          "X10193hl_28_HC5.2"         
[89] "X10193hl_3_NC1.1"           "X10193hl_35_NC2.3"         
[91] "X10193hl_38_NC3.3"          "X10193hl_39_HC3.3"         
[93] "X10193hl_40_NC4.3"          "X10193hl_42_NC5.3"         
[95] "X10193hl_43_HC5.3"          "X10193hl_5_NC2.1"          
[97] "X10193hl_6_HC2.1"           "X10193hl_7_NC3.1"          
[99] "X10193hl_8_HC3.1"           "X10193hl_9_NC4.1"          
[101] "X10194hl_10_HC4.4"          "X10194hl_11_NC5.4"         
[103] "X10194hl_12_HC5.4"          "X10194hl_14_NC6.4"         
[105] "X10194hl_18_NC1.5"          "X10194hl_20_NC2.5"         
[107] "X10194hl_22_NC3.5"          "X10194hl_23_HC3.5"         
[109] "X10194hl_24_NC4.5"          "X10194hl_25_HC4.5"         
[111] "X10194hl_27_NC5.5"          "X10194hl_28_HC5.5"         
[113] "X10194hl_3_NC1.4"           "X10194hl_35_NC2.6"         
[115] "X10194hl_36_HC2.6"          "X10194hl_38_NC3.6"         
[117] "X10194hl_39_HC3.6"          "X10194hl_40_NC4.6"         
[119] "X10194hl_41_HC4.6"          "X10194hl_42_NC5.6"         
[121] "X10194hl_43_HC5.6"          "X10194hl_5_NC2.4"          
[123] "X10194hl_6_HC2.4"           "X10194hl_7_NC3.4"          
[125] "X10194hl_8_HC3.4"           "X10194hl_9_NC4.4"          
> 
  > Tr.clu2 <- Tr.PAM[Tr.PAM[,ncol(Tr.PAM)]==2,]
> rownames(Tr.clu2)
[1] "X13122mo_62_WT_A_P1_24h_LC" "X10201hl_15_NC7.1"          "X10201hl_16_VC7.1"         
[4] "X10201hl_30_NC7.2"          "X10201hl_31_VC7.2"          "X10201hl_45_NC7.3"         
[7] "X10201hl_46_VC7.3"          "X10202hl_15_NC7.4"          "X10202hl_16_VC7.4"         
[10] "X10202hl_30_NC7.5"          "X10202hl_31_VC7.5"          "X10202hl_44_VC6.6"         
[13] "X10202hl_45_NC7.6"          "X10202hl_46_VC7.6"          "X10193hl_15_HC6.1"         
[16] "X10193hl_16_NC7.1"          "X10193hl_17_HC7.1"          "X10193hl_31_NC7.2"         
[19] "X10193hl_32_HC7.2"          "X10193hl_45_HC6.3"          "X10193hl_46_HC7.3"         
[22] "X10193hl_47_NC7.3"          "X10194hl_16_NC7.4"          "X10194hl_17_HC7.4"         
[25] "X10194hl_31_NC7.5"          "X10194hl_32_HC7.5"          "X10194hl_46_HC7.6"         
[28] "X10194hl_47_NC7.6" 

### Average silhouette coefficient to choose number of clusters for PAM of RV matrix ###

> num.clus(rv_matrix)
2  avg:  0.13 clus:  0.1 0.26 
3  avg:  0.17 clus:  0.09 0.19 0.4 
4  avg:  0.16 clus:  0.06 0.25 0.13 0.34 
5  avg:  0.17 clus:  0.05 0.23 0.15 0.15 0.37 
6  avg:  0.18 clus:  0.04 0.26 0.45 0.12 0.15 0.37 
7  avg:  0.22 clus:  0.95 0.25 0.45 0.12 0.05 0.14 0.37 
8  avg:  0.23 clus:  0.95 0.26 0.25 0.09 0.1 0.29 0.13 0.44 
9  avg:  0.24 clus:  0.95 0.26 0.44 0.08 0.04 0.07 0.44 0.36 0.2 
10  avg:  0.25 clus:  0.95 0.26 0.42 0.11 0 0.08 0.22 0.15 0.44 0.32 
11  avg:  0.24 clus:  0.95 0.26 0.42 0.22 -0.01 0.11 0.22 0.07 0.14 0.43 0.45 
12  avg:  0.25 clus:  0.95 0.26 0.42 0.21 -0.01 -0.03 0.22 0.22 0.43 0.39 0.18 0.25 
13  avg:  0.28 clus:  0.95 0.26 0.84 0.06 -0.03 0.03 0.22 0.22 0.39 0.41 0.39 0.17 0.25 
14  avg:  0.28 clus:  0.95 0.26 0.84 0.06 -0.03 0.08 0.15 0.21 0.37 0.18 0.41 0.39 0.21 0.25 
15  avg:  0.28 clus:  0.95 0.26 0.84 0.05 0 0.07 0.15 0.08 0.19 0.37 0.18 0.41 0.37 0.21 0.25 
16  avg:  0.27 clus:  0.95 0.25 0.84 0.05 0 0.07 0.15 0.08 0.17 0.36 0.17 0.41 0.37 0.22 0 0.25 
17  avg:  0.26 clus:  0.95 0.25 0.84 0.05 0 0.06 0.15 0 0.17 0.33 0.17 0.41 0 0.37 0.22 0 0.23 
18  avg:  0.27 clus:  0.95 0.25 0.84 0.05 0 0.2 0 0.15 0 0.17 0.33 0.17 0.41 0 0.37 0.22 0 0.23 
19  avg:  0.25 clus:  0.95 0.22 0.84 0.05 0 0.2 0 0 0 0.17 0.33 0.14 0.41 0 0.37 0.22 0 0.23 0 
20  avg:  0.24 clus:  0.95 0.22 0.84 0.05 0 0.2 0 0 0 0.12 0.33 0.14 0.41 0 0.34 0.22 0 0.33 0 0

### Samples within PAM clusters of RV matrix ###

#10 clusters according to silhouette coefficient
> clu1 <- PAM[PAM[,ncol(PAM)]==1,]
> rownames(clu1)
[1] "Syn exp1 HC" "Syn exp3 HC"

> clu2 <- PAM[PAM[,ncol(PAM)]==2,]
> rownames(clu2)
[1] "Syn exp1 3h LC"      "Syn exp3 3h LC"      "Ara exp4 light HC"   "Ara exp5 dark 1d LC"

> clu3 <- PAM[PAM[,ncol(PAM)]==3,]
> rownames(clu3)
[1] "Syn exp1 24h LC" "Syn exp2 HC"     "Syn exp3 24h LC"

> clu4 <- PAM[PAM[,ncol(PAM)]==4,]
> rownames(clu4)
[1] "Syn exp2 3h LC"  "Ara exp6 6h NC"  "Ara exp6 8h NC"  "Ara exp6 10h NC"
[5] "Ara exp6 2h VC" 

> clu5 <- PAM[PAM[,ncol(PAM)]==5,]
> rownames(clu5)
[1] "Syn exp2 24h LC"      "Ara exp4 light 3d LC" "Ara exp4 light 5d LC"
[4] "Ara exp6 12h NC"      "Ara exp6 4h VC"      

> clu6 <- PAM[PAM[,ncol(PAM)]==6,]
> rownames(clu6)
[1] "Ara exp4 light 1d LC" "Ara exp5 3d LC"       "Ara exp7 2h NC"      
[4] "Ara exp7 2h HC"      

> clu7 <- PAM[PAM[,ncol(PAM)]==7,]
> rownames(clu7)
[1] "Ara exp5 dark HC"    "Ara exp5 dark 5d LC" "Ara exp6 2h NC"      "Ara exp6 6h VC"     
[5] "Ara exp7 4h HC"     

> clu8 <- PAM[PAM[,ncol(PAM)]==8,]
> rownames(clu8)
[1] "Ara exp6 4h NC"  "Ara exp6 12h VC" "Ara exp7 4h NC"  "Ara exp7 6h NC" 
[5] "Ara exp7 12h NC" "Ara exp7 12h HC"

> clu9 <- PAM[PAM[,ncol(PAM)]==9,]
> rownames(clu9)
[1] "Ara exp6 24h NC" "Ara exp6 24h VC" "Ara exp7 24h NC" "Ara exp7 24h HC"

> clu10 <- PAM[PAM[,ncol(PAM)]==10,]
> rownames(clu10)
[1] "Ara exp6 8h VC"  "Ara exp6 10h VC" "Ara exp7 8h NC"  "Ara exp7 10h NC"
[5] "Ara exp7 6h HC"  "Ara exp7 8h HC"  "Ara exp7 10h HC"

#4 clusters according to pearson clustering
> clu1 <- PAM[PAM[,ncol(PAM)]==1,]
> rownames(clu1)
[1] "Syn exp1 HC"          "Syn exp2 HC"          "Syn exp2 3h LC"      
[4] "Syn exp2 24h LC"      "Syn exp3 HC"          "Ara exp4 light 1d LC"
[7] "Ara exp4 light 3d LC" "Ara exp4 light 5d LC" "Ara exp5 dark HC"    
[10] "Ara exp5 dark 1d LC"  "Ara exp5 3d LC"       "Ara exp5 dark 5d LC" 
[13] "Ara exp6 2h NC"       "Ara exp6 12h NC"      "Ara exp6 4h VC"      
[16] "Ara exp6 6h VC"       "Ara exp7 2h HC"       "Ara exp7 4h HC"      

> clu2 <- PAM[PAM[,ncol(PAM)]==2,]
> rownames(clu2)
[1] "Syn exp1 3h LC"    "Syn exp1 24h LC"   "Syn exp3 3h LC"    "Syn exp3 24h LC"  
[5] "Ara exp4 light HC" "Ara exp6 8h NC"    "Ara exp6 24h NC"   "Ara exp6 24h VC"  
[9] "Ara exp7 24h NC"   "Ara exp7 24h HC"  

> clu3 <- PAM[PAM[,ncol(PAM)]==3,]
> rownames(clu3)
[1] "Ara exp6 4h NC"  "Ara exp6 6h NC"  "Ara exp6 10h NC" "Ara exp6 12h VC"
[5] "Ara exp7 2h NC"  "Ara exp7 4h NC"  "Ara exp7 6h NC"  "Ara exp7 12h NC"
[9] "Ara exp7 12h HC"

> clu4 <- PAM[PAM[,ncol(PAM)]==4,]
> rownames(clu4)
[1] "Ara exp6 2h VC"  "Ara exp6 8h VC"  "Ara exp6 10h VC" "Ara exp7 8h NC" 
[5] "Ara exp7 10h NC" "Ara exp7 6h HC"  "Ara exp7 8h HC"  "Ara exp7 10h HC"

### Average silhouette coefficient to choose number of clusters for PAM of RV matrix ###
### No "dark" Arabidopsis samples in dataset ###

> num.clus(rv_matrix1)
2  avg:  0.19 clus:  0.12 0.41 
3  avg:  0.16 clus:  0.09 0.09 0.4 
4  avg:  0.19 clus:  0.06 0.19 0.22 0.38 
5  avg:  0.22 clus:  0.38 0.18 0.2 0.1 0.37 
6  avg:  0.26 clus:  0.95 0.16 0.31 0.18 0.1 0.37 
7  avg:  0.27 clus:  0.95 0.41 0.47 0.18 0.01 0.07 0.37 
8  avg:  0.27 clus:  0.95 0.41 0.44 0.01 0.11 0.13 0.19 0.34 
9  avg:  0.26 clus:  0.95 0.41 0.44 0.01 0.11 0.13 0.19 0.38 0.24 
10  avg:  0.26 clus:  0.95 0.41 0.45 0.22 -0.03 0.09 0.24 0.09 0.38 0.24

