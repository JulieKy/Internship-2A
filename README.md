# Preterm Neonates Breathing Sounds Analysis 

## Abstract 
Combining signal processing and biomedical, the project focuses on the breathing sound analysis of premature newborns with respiratory deficiency. In connection with a medical team, the work consisted first in signal preprocessing to remove the noises (crying, heartbeat, machine noises, ...) by relying on labeled signals. Then, the extraction of signal characteristics was necessary to learn more about Surfactant Replacement Therapy. 

## Vocabulary
**CS:** Crying sections 

**NCS:** Non-crying sections

**RDS:** Respiratory Distress Syndrome

**SRT:** Surfactant Replacement Therapy

## File Organization

### CODES 
``` main_Julie.m``` : Main file of the project. Deals with functions to do the preprocessing and fetaures extraction, and creates and Excel file containing the different features.

#### Preprocessing
##### Crying removing 
``` crying_learning.m``` : Main file for learning the features of CS and NCS.

``` labelling.m``` : Labels the recordings using Audacity text files, thank to a moving average amd give the Fleiss'es KAPPA coefficients.

``` Fleiss.m``` : ompute the Fleiss'es kappa.

``` label2signal.m```: With the placement of the CS/NCS in the label_final matrix returned by LABELLING.m, this function gives the corresponding piece of signal.

``` power_ratio_band.m```: Calculate the power ratio of NCS/CS (depending on flag_section) of a signal. 

``` spectrogram_CS.m```: Display the spectrogram and the time representation with CS and NCS of a part of a signal AND return its spectrogram. 

``` NCS_CS_features_boxplot.m```: Extract some features of CS and NCS, return them and do their boxplots.

``` crying_removing.m```: Remove the CS in all the signals
##### Filtering 
```filterbp.m ```: Filter signals using the ButterWorth Filter. 


#### Features 
```spectral_features.m```: Main file of feature extraction. Computes the following spectral features: *meanPSD; stdPSD; medPSD; bw; p25; p75; IQR; TP; p100_200; p200_400; p400_800; spectrum_slope2; r_square2; nb_pks_MAF;  f_higherPk_MAF; dif_higherPks_MAF; nb_pks_GMM;  f_higherPk_GMM; dif_higherPks_GMM; GMM_parameters; pxx*.

```Welch_periodogram.m```: Computes the Welch Periodogram.

```Welch_periodogram.m```: Computes peaks features: *higher peak; number of peaks; their frequencies*.

```lpc_lsf_coeff.m```: Gives the first 6 coefficients of LPC and LSF.

```temporal_features.m```:Computes the Zero Crossing Rate.

```AUTRES_FEATURES.m```: A FAIRE!!


#### Display
```display_NCS_CS_annotations.m```: Displays the annotated labels CS and NCS of a signal. 

```display_PR_NCS_CS_interquartiles.m```: Displays the periodograms of annotated NCS and CS. 

### DATA

#### Labels 
Audacity text files corresponding of the CS/NCS annotations of three observators. There are named like this: ```ObservatorID_SampleID```.

#### Samples
Samples used in a previous project. 

#### Samples_Belle 
37 samples of preterm neonates with RDS. They are used in this project. The recordings have been done on babies who need SRT, before and after taking surfactant, exactly on Day 2 and on Day 28 after birth. 

### FIGURES

#### Preprocessing figures 
Gathers all the preprocessing figures (labelling, crying detection, power ratio of CS/NCS, spectrograms of CS/NCS, boxplot of CS/NCS features, ...)

#### Time figures 
Gathers the time representation of all samples. 