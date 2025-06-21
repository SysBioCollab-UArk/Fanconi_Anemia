from fanconi_anemia import model
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import *
import os

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

time_perturb_value = [{0: ('Lesion(fanc=None, ner=None)', 44 * 6.4)},  # Expt A
                      {0: ('ICL(b=None)', 2 * 6.4)},  # Expt B
                      {0: [('Lesion(fanc=None, ner=None)', 5.3 * 6.4), ('ICL(b=None)', 38.7 * 6.4)]}]  # Expt C

sim_protocols = [SequentialInjections(solver, t_equil=100, time_perturb_value=tpv) for tpv in time_perturb_value]

custom_priors = {}  # {'N': ('uniform', 0.3)}
no_sample = ['ICL_0', 'DSB_0', 'Lesion_0', 'k_ICL_synth', 'k_AG20_lump', 'k_BL100_lump', 'k_CEF_lump',
             'k_AG20_BL100_CEF_lump']
obs_labels = {'DNA_lesions': 'MAs', 'Interstrand_crosslinks': 'ICLs'}

exp_data_file = os.path.join('DATA', 'Averbeck1988_Normal.csv')

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocols,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': True, 'save_sim_data': True})
