Matlab code for Negroni et al. rabbit ventricular model with myofilament contraction.

This model describes excitation-contraction coupling in the rabbit ventricular myocytes.
A new model of myofilament contraction is here included into an established (and updated) 
computational framework, which integrates descriptions of electrophysiology, Ca and Na
handling (from Shannon et al. Biophys J. 2004 Nov;87(5):3351-71), Ca/calmodulin-dependent
protein kinase II (CaMKII) and protein kinase A (PKA) signaling pathways (from Soltis and
Saucerman, Biophys J. 2010 Oct 6;99(7):2038-47).

______________________________________________________________________________________
Contents:

readme.txt				this file

.m files

rabbit_myofilament_masterCompute.m	loads initial conditions and runs the simulation
rabbit_myofilament_masterODEfile.m	integrates the following model components
rabbit_myofilament_barODEfile.m		beta-adrenergic (PKA) phosphorylation module
rabbit_myofilament_camkiiODEfile.m	CaMKII phosphorylation module
rabbit_myofilament_camODEfile.m		CaM module
rabbit_myofilament_eccODEfile.m		excitation-contraction coupling module

.mat files				initial conditions (obtained at 1 Hz pacing)

- yfinal_myo_WT_1Hz_isometric		WT model (isometric contraction)
- yfinal_myo_WT_1Hz_isometric_240iso	WT model + 240 s ISO administration (0.1 uM)
- yfinal_myo_WT_1Hz_isometric_240iso_X	WT model + 240 s ISO administration (0.1 uM)
					with reduced PKA effect on target X

- yfinal_myo_WT_1Hz_isotonic		WT model (isotonic contraction)
- yfinal_myo_WT_1Hz_isotonic_240iso	WT model + 240 s ISO administration (0.1 uM)
- yfinal_myo_WT_1Hz_isotonic_240iso_X	WT model + 240 s ISO administration (0.1 uM)
					with reduced PKA effect on target X
______________________________________________________________________________________


Reference:

J.A. Negroni, S. Morotti, E.C. Lascano, A.V. Gomes, E. Grandi, J.L. Puglisi, D.M. Bers.
ß-adrenergic effects on cardiac myofilaments and contraction in an integrated rabbit
ventricular myocyte model.
J Mol Cell Cardiol. 2015 Apr; 81:162-75.
doi: https://doi.org/10.1016/j.yjmcc.2015.02.014

Please cite the above paper when using this model.
