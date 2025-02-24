from fanconi_anemia import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *
import os

# TODO: Comment out unnecessary Observables

print(model.observables)
quit()

# Monomer('FANCA', ['fancg', 'faap20']),
# Monomer('FANCG', ['fanca']),
# Monomer('FAAP20', ['fanca']),
# Monomer('FANCB', ['fancl', 'faap100']),
# Monomer('FANCL', ['fancb', 'faap100']),
# Monomer('FAAP100', ['fancb', 'fancl']),
# Monomer('FANCF', ['fancc']),
# Monomer('FANCC', ['fancf', 'fance']),
# Monomer('FANCE', ['fancc']),
# Monomer('FANCM', ['dna', 'facpx']),
# Monomer('FANCT', ['facpx']),
# Monomer('FANCI', ['fancd2', 'fancp', 'state'], {'state': ['x', 'ub']}),
# Monomer('FANCD2', ['fanci', 'facpx', 'fancp', 'dna', 'state'], {'state': ['x', 'ub']}),
# Monomer('FANCP', ['fanci', 'fancd2', 'fancq']),
# Monomer('FANCQ', ['fancp']),

# Below are reasonable starting guesses for initial concentrations of key FA pathway proteins, assuming a
# nuclear volume of ~300 fL for a human epithelial cell:
# -------------------------------------------------------
# Protein	        Estimated Copy Number	Estimated Concentration (nM)	Source/Justification
# FANCA	            ~200,000	            ~1,000 nM (1 μM)	            Abundant core complex component
# FANCB	            ~50,000	                ~250 nM	                        Lower expression than FANCA
# FANCC	            ~100,000	            ~500 nM	                        Core complex component
# FANCE	            ~50,000	                ~250 nM	                        Facilitates FANCD2 activation
# FANCF	            ~50,000	                ~250 nM	                        Core component
# FANCG	            ~75,000	                ~375 nM	                        Interacts with FANCA
# FANCL	            ~10,000 - 25,000	    ~50 - 125 nM	                E3 ubiquitin ligase, lower abundance
# FAAP20	        ~10,000 - 20,000	    ~50 - 100 nM	                Stabilizes FANCA
# FAAP100	        ~20,000 - 30,000	    ~100 - 150 nM	                Core complex component
# FANCM	            ~50,000 - 100,000	    ~250 - 500 nM	                DNA translocase, binds stalled forks
# FANCD2	        ~100,000 - 200,000	    ~500 - 1,000 nM	                Central mediator of repair
# FANCI	            ~100,000 - 200,000	    ~500 - 1,000 nM	                Forms heterodimer with FANCD2
# FANCT (UBE2T)	    ~10,000 - 50,000	    ~50 - 250 nM	                E2 ubiquitin-conjugating enzyme
# FANCP (SLX4)	    ~10,000 - 50,000	    ~50 - 250 nM	                Endonuclease scaffold
# FANCQ (XPF/ERCC4)	~50,000 - 100,000	    ~250 - 500 nM	                Nuclease involved in unhooking
# ERCC1             ~18,000 - 90,000        ~10 - 50 nM

# Experiment A: 44 * 6.4 = 281.6 DNA adducts, 0 ICLs
# Experiment B: 0 DNA adducts, 2 * 6.4 = 12.8 ICLs
# Experiment C: 5.3 * 6.4 = 33.92 DNA adducts, 38.7 * 6.4 = 247.68 ICLs

solver = ScipyOdeSimulator(model)

time_perturb_value = [{0: ('Lesion(fanc=None, ner=None)', 44 * 6.4)},  # Expt A, Normal
                      {0: ('ICL(b=None)', 2 * 6.4)},  # Expt B, Normal
                      {0: [('Lesion(fanc=None, ner=None)', 5.3 * 6.4), ('ICL(b=None)', 38.7 * 6.4)]},  # Expt C, Normal
                      {0: ('Lesion(fanc=None, ner=None)', 44 * 6.4)},  # Expt A, FA 150
                      {0: ('ICL(b=None)', 2 * 6.4)},  # Expt B, FA 150
                      {0: [('Lesion(fanc=None, ner=None)', 5.3 * 6.4), ('ICL(b=None)', 38.7 * 6.4)]},] # Expt C, FA 150

multi_exp_injection = ParallelExperiments(solver, t_equil=1e3, time_perturb_value=time_perturb_value)

custom_priors = {}  # {'N': ('uniform', 0.3)}
no_sample = ['ICL_0', 'DSB_0', 'Lesion_0', 'k_ICL_synth', 'k_AG20_lump', 'k_BL100_lump', 'k_CEF_lump',
             'k_AG20_BL100_CEF_lump']
obs_labels = {'DNA_lesions': 'MAs', 'Interstrand_crosslinks': 'ICLs'}

exp_data_file = os.path.join('DATA', 'Averbeck1988.csv')

params = ['FANCA_0', 'kf_AG', 'kr_AG', 'kf_A20', 'kr_A20', 'kf_AG20_BL100', 'kr_AG20_BL100', 'kf_AG20_BL100_CEF',
 'kr_AG20_BL100_CEF', 'kf_AG20_CEF', 'kr_AG20_CEF', 'kf_AG20_CEF_BL100', 'kr_AG20_CEF_BL100', 'kf_FAcpx_M',
 'kr_FAcpx_M', 'kf_FAcpx_T', 'kr_FAcpx_T', 'kf_FAcpxMT_binds_ID2', 'kr_FAcpxMT_binds_ID2',
 'k_FAcpxMT_ICL_release_ID2ub', 'k_FAcpxMT_Lesion_release_ID2ub']

param_expts_map = {param: [('A', 'B', 'C'), ('AA', 'BB', 'CC')] for param in params}

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      multi_exp_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample,
                                      param_expts_map=param_expts_map)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': True, 'save_sim_data': True})
