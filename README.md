# AssociativeMotorAdaptation
## Codes and data accompanying the paper: Avraham et al. Contextual effects in sensorimotor adaptation adhere to associative learning rules.
eLife 2022;11:e75801. DOI: https://doi.org/10.7554/eLife.75801

### Codes:

  DataAnalysis:
  This folder includes custom-written Matlab codes (.m) to analyze the data from Experiments 1-4. The name of each code file ends with "main".
  The folder also contains Matlab data files (.mat) that can be generated by the data analysis codes.
 
  Simulations:
  This folder includes simulation codes (.m) and simulated data (.mat)
 
  Functions:
  Contains all custom-written functions for statistical analyses and plotting (called by data analysis and simulation codes).

### Data:
  
  Includes raw data collected during the experiments:
    
    Experiment 1- AssociativeAdaptation_Exp1_Differential_trials.mat
    Experiment 2- AssociativeAdaptation_Exp2_DifferentialTiming_Delay_trials.mat
    Experiment 3- AssociativeAdaptation_Exp3_DifferentialTiming_Control_trials.mat
    Experiment 4- AssociativeAdaptation_Exp4_Compound_trials.mat

    Each data file contain a Table T with the following variables (variables can change across experiments):
    SN- subject number
    tester- experimenter initial
    group
    cond- whether frame (f, light) or tone (t) was associated with CS+
    TN- trial number
    CN- cycle number
    BN- block number
    CCW- rotation direction during learning: 1- Counterclockwise; 0- Clockwise
    tgtsize- target size (mm diameter)
    hand_theta- hand movement direction with respect to the target direction at the radial distance to the target (this variable was used as the 'hand angle' measure in the paper)
    hand_theta_maxv- hand movement direction with respect to the target direction at maximum movement velocity
    hand_theta_maxradv- hand movement direction with respect to the target direction at maximum movement radial velocity
    handMaxRadExt- hand movement direction with respect to the target direction at maximum radial extent
    hand_theta_50- hand movement direction with respect to the target direction at 50 msec after movement initiation
    raw_ep_hand_ang- actual hand movement direction at the radial distance to the target (without subtracting target location)
    ti- target location (deg)
    fbi- feedback: 1- cursor presented; 0- no cursor
    ri- rotation size and direction: positive values- counterclockwise rotation; negative values- clockwise rotation.
    clampi- clamped feedback: 1- the rotation is with respect to the target location (Experiments 2 and 2S); 0- the rotation is with respect to the hand (Experiment 1)
    cs_tone/tone- CS tone: 1- present; 0- absent 
    cs_frame/light- CS light: 1- present; 0- absent
    cs_time- interval between trial onset and cs onset (s)
    go_time- interval between trial onset and imperative (s)
    MT- movement time (the interval between the time at which the amplitude of the movement exceeded 1 cm from the start location to the time at which the amplitude reached the radial distance of the target)
    RT- reaction time (the interval between the imperative onset and the time that the hand position exceeded a distance of 1 cm from the start location)
    ST- search time (the time required to move back to the start location)
    radvelmax- maximum radial velocity
    maxRadDist- maximum radial movement extent

    More details on how the data was processed and analyzed can be found in the Methods section of the paper.
