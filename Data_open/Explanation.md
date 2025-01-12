Explanation of Data in this folder
data_bkg1:
'file':file name
'frame_size':size of frame (30×20)
'stim': position of red dot (Trials×12(dot repetition)×2(X,Y))
'frame':position of background rectanele (Trials×2(X,Y))
'fixation':position of fixation cross (Trials×2(X,Y))
'fp_on':timing of fixation cross presentation relative to XXX.
'xyz':position of recording site (LR, AP, DV)
'dep':depth of recording site from surface 
'ap':anterior-posterior position relative to center of grid which was used at recording
'angle':angle of penetration
'monkey':identification of monkey monkey1Right:4 monkey1Left4.5 Monkey2Left:5 Monkey3Left:2
'spike_ts':spike timistamp
'photo_ts':timing of dot presentation detected by photo sensor attached monitor


data_nobkg:
data with background (data_bkg_0_1) or withoug background (data_bkg_0_0).
this is used for Figure 8.

Entropy_data.mat:
(number(50 or 100) is the time window(msec))
Entropy:neural informtion 
I_SHUFF:neural information which was calculated by shuffled background location across trials
P_CHIV:P valuue which was culated by shuffled background location (H2)
P_CHIV_dotsh:P valuue which was culated by shuffled dot location (H1)
spike_dot_bkg_filter:P valuue which was culated by shuffled trials
bkgpos_map: numnber of background presentation aligned with retinotopic coordinate
spike_peridot_full:number of spikes aligned with background position 

Entropy_data_nobkg:entropy data without background presentation

meanS and meanSNo: average number of spikes for each dot presentation 

anatomy.mat
Position and brain are of each neuron 
brain area 1:V2 2:V6A 3:PEc 4:7M 5:MIP

cell_POS.mat
cell_XY_all anatomical location of cells
cell_POS_all anatomical location of PO sulcus


gaze_data
position of gaze (167973Trials: Number of Trials×Number of Neuron)
psth:PeriStimulusSpikeCounts of each trial, each Dot presentation, sorted fixation presentation


gabor_bkg_abs_each_0005:
result of gabor fitting and fitting parameters for dot/bkg neurons (gabor_bkg_abs_each_0005)
and dot ret neurons (gabor_ret_abs_each_0005)

hybrid_all
Entropy_thre:neural information
R_co:Correlation coefficients representing dot/ret (ordinate) and dot/bkg (abscissa) information in the
hybrid coding analysis


trans_bkg
GchiV:The chi-square values calculated separately when the the bakcground frame appeared in quadrants 1-4
Spike_dot_filter_g:receptive fields map calculated separately when the the bakcground frame appeared in quadrants 1-4


gain_bkg(retino)
GchiV :The chi-square values calculated separately when the fixation point appeared in quadrants 1-4


Entropy_data_binstart_0dotbkg_100(50)
Entrropy:neural information evoked by each dot presentation (dot1~12)
I_SHUFF:neural information from shuffled data evoked by each dot presentation (dot1~12)
