/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/bio_opt.h
 *  
 *  Description:
 *  Header file for bio_opt
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: bio_opt.h 7195 2022-09-13 10:34:31Z bai155 $
 *
 */

#if !defined(_BIO_OPT_H)
#define _BIO_OPT_H

/*
 * Bio-optical properties
 */
typedef struct {
  double* landa;
  double* kw_s;
  double* kw_s2;
  double* bw_s;
  double* bw_s2;
  double* ay_A;
  double* ay_A2;
  double* ay_B;
  double* ay_B2;
  double* ad_A;
  double* ad_A2;
  double* ad_B;
  double* ad_B2;
  double* bbp_A;
  double* bbp_A2;
  double* bbp_B;
  double* bbp_B2;
  double* ad_CarbSand_A;
  double* ad_CarbSand_A2;
  double* ad_CarbSand_B;
  double* ad_CarbSand_B2;
  double* bbp_CarbSand_A;
  double* bbp_CarbSand_A2;
  double* bbp_CarbSand_B;
  double* bbp_CarbSand_B2;
  double* scatfrac;
  double* MA_aAwave;double* MAR_aAwave;double* MAG_aAwave;
  double* SG_aAwave;
  double* SGH_aAwave;
  double* MA_aAwave2;double* MAR_aAwave2;double* MAG_aAwave2;
  double* SG_aAwave2;
  double* SGH_aAwave2;
  double gone;
  double gtwo;
  double acdom443star;
  double *sand_srs;
  double *mud_srs;
  double *finesed_srs;
  double *coral_skeleton_srs;
  double *macroalgae_srs;
  double *seagrass_srs;
  double *yC_diadinoxanthin;
  double *yC_diatoxanthin;
  double *yC_chlorophylla;
  double *yC_chlorophyllc2;
  double *yC_peridinin;
  double *yC_betacarotene;
  double *yC_zeaxanthin;
  double *yC_chlorophyllb;
//    double *yC_divinyl_chlb;
  double *yC_alloxanthin;
  double *yC_19_BF;
  double *yC_19_HF;
  double *yC_myxoxanthophyll;
  double *yC_fucoxanthin;
  double *yC_echinenone;
  double *yC_PE;
  double *yC_PC;
  double *yC_allophycocyanin;
  double *yC_picoplankton;
  double *yC_rhodomonas_duplex;
  double *yC_microplankton;
  double *yC_tricho;
  double *yC_symbiodinium;
  double *yC_diadinoxanthin_rsr;
  double *yC_diatoxanthin_rsr;
  double *yC_chlorophylla_rsr;
  double *yC_chlorophyllc2_rsr;
  double *yC_peridinin_rsr;
  double *yC_betacarotene_rsr;
  double *yC_zeaxanthin_rsr;
  double *yC_chlorophyllb_rsr;
//    double *yC_divinyl_chlb_rsr;
  double *yC_myxoxanthophyll_rsr;
  double *yC_alloxanthin_rsr;
  double *yC_19_BF_rsr;
  double *yC_19_HF_rsr;
  double *yC_fucoxanthin_rsr;
  double *yC_echinenone_rsr;
  double *yC_PE_rsr;
  double *yC_PC_rsr;
  double *yC_allophycocyanin_rsr;
  double *yC_picoplankton_rsr;
  double *yC_rhodomonas_duplex_rsr;
  double *yC_microplankton_rsr;
  double *yC_tricho_rsr;
  double *yC_symbiodinium_rsr;
  double *aC_AUS1;
  double *aC_AUS2;
  double *aC_ICE1;
  double *aC_ICE2;
  double *aC_ICE3;
  double *aC_KUW1;
  double *aC_KUW2;
  double *aC_NIG1;
  double *aC_SAH1;
  double *aC_SAH2;
  double *aC_OAH1;
  double *aC_OAH2;
  double *aC_CAL1;
  double *aC_CAL2;
  double *aC_QUA1;
  double *aC_ILL1;
  double *aC_ILL2;
  double *aC_KAO1;
  double *aC_KAO2;
  double *aC_KAO3;
  double *aC_MON1;
  double *aC_MON2;
  double *aC_SAN1;
  double *bC_AUS1;
  double *bC_AUS2;
  double *bC_ICE1;
  double *bC_ICE2;
  double *bC_ICE3;
  double *bC_KUW1;
  double *bC_KUW2;
  double *bC_NIG1;
  double *bC_SAH1;
  double *bC_SAH2;
  double *bC_OAH1;
  double *bC_OAH2;
  double *bC_CAL1;
  double *bC_CAL2;
  double *bC_QUA1;
  double *bC_ILL1;
  double *bC_ILL2;
  double *bC_KAO1;
  double *bC_KAO2;
  double *bC_KAO3;
  double *bC_MON1;
  double *bC_MON2;
  double *bC_SAN1;
  double *LJCO_aNAP;
  double *LJCO_bNAP;
  double *LJCOc_aNAP;
  double *LJCOc_bNAP;
  double *aC_AUS1_rsr;
  double *aC_AUS2_rsr;
  double *aC_ICE1_rsr;
  double *aC_ICE2_rsr;
  double *aC_ICE3_rsr;
  double *aC_KUW1_rsr;
  double *aC_KUW2_rsr;
  double *aC_NIG1_rsr;
  double *aC_SAH1_rsr;
  double *aC_SAH2_rsr;
  double *aC_OAH1_rsr;
  double *aC_OAH2_rsr;
  double *aC_CAL1_rsr;
  double *aC_CAL2_rsr;
  double *aC_QUA1_rsr;
  double *aC_ILL1_rsr;
  double *aC_ILL2_rsr;
  double *aC_KAO1_rsr;
  double *aC_KAO2_rsr;
  double *aC_KAO3_rsr;
  double *aC_MON1_rsr;
  double *aC_MON2_rsr;
  double *aC_SAN1_rsr;
  double *bC_AUS1_rsr;
  double *bC_AUS2_rsr;
  double *bC_ICE1_rsr;
  double *bC_ICE2_rsr;
  double *bC_ICE3_rsr;
  double *bC_KUW1_rsr;
  double *bC_KUW2_rsr;
  double *bC_NIG1_rsr;
  double *bC_SAH1_rsr;
  double *bC_SAH2_rsr;
  double *bC_OAH1_rsr;
  double *bC_OAH2_rsr;
  double *bC_CAL1_rsr;
  double *bC_CAL2_rsr;
  double *bC_QUA1_rsr;
  double *bC_ILL1_rsr;
  double *bC_ILL2_rsr;
  double *bC_KAO1_rsr;
  double *bC_KAO2_rsr;
  double *bC_KAO3_rsr;
  double *bC_MON1_rsr;
  double *bC_MON2_rsr;
  double *bC_SAN1_rsr;
  double *LJCO_aNAP_rsr;
  double *LJCO_bNAP_rsr;
  double *LJCOc_aNAP_rsr;
  double *LJCOc_bNAP_rsr;

  double *B_carb_rsr;
  double *B_terr_rsr;

  double *a0;
  double *a1;
  double *a2;
  double *a3;
  double *b1;
  double *b2;
  double *b3;
  double *d1;
  double *d2;
  double *d3;


} bio_opt_prop;


bio_opt_prop *bio_opt_init(ecology* e);
void bio_opt_free(bio_opt_prop *b);

#endif
