from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# STEPS:
# 1. Parp1 binds DSB
# 2. MRE11 binds Parp1
# 3. CtIP binds Parp1


Monomer("DSB", ["b"])
Monomer("Parp1", ["dsb", "ctip_mre11"])
Monomer("CtIP", ["parp1"])
Monomer("MRE11", ["parp1"])

Parameter("DSB_0", 100)
Parameter("Parp1_0", 100)
Parameter("CtIP_0", 100)
Parameter("MRE11_0", 100)

Initial(DSB(b=None), DSB_0)
Initial(Parp1(dsb=None, ctip_mre11=None), Parp1_0)
Initial(CtIP(parp1=None), CtIP_0)
Initial(MRE11(parp1=None), MRE11_0)

Parameter("kf_DSB_Parp1", 1)
Parameter("kr_DSB_Parp1", 10)
Parameter("kf_DSB_CtIP", 1)
Parameter("kr_DSB_CtIP", 10)

# Monomer definitions for two distinct Parp1 molecules
Monomer("Parp1_CtIP", ["dsb", "ctip"])
Monomer("Parp1_MRE11", ["dsb", "mre11"])

# Rule for MRE11 binding first to Parp1_MRE11
Rule("Parp1_MRE11_binds_DSB",
     DSB(b=None) + Parp1_MRE11(dsb=None, mre11=None) |
     DSB(b=1) % Parp1_MRE11(dsb=1, mre11=None), kf_DSB_Parp1, kr_DSB_Parp1)

Rule("MRE11_binds_Parp1",
     MRE11(parp1=None) + DSB(b=1) % Parp1_MRE11(dsb=1, mre11=None) |
     MRE11(parp1=2) % DSB(b=1) % Parp1_MRE11(dsb=1, mre11=2), kf_DSB_Parp1, kr_DSB_Parp1)

# Rule for CtIP binding later to Parp1_CtIP
Rule("Parp1_CtIP_binds_DSB",
     DSB(b=None) + Parp1_CtIP(dsb=None, ctip=None) |
     DSB(b=1) % Parp1_CtIP(dsb=1, ctip=None), kf_DSB_Parp1, kr_DSB_Parp1)

Rule("CtIP_binds_Parp1",
     CtIP(parp1=None) + DSB(b=1) % Parp1_CtIP(dsb=1, ctip=None) |
     CtIP(parp1=2) % DSB(b=1) % Parp1_CtIP(dsb=1, ctip=2), kf_DSB_CtIP, kr_DSB_CtIP)

# Next step in which the resection process begins to create ssDNA.
class ResectionProcess:
     def __init__(self, ssDNA_overhang_length, ctip_active, mre11_active, exo1_active, dna2_active):
          self.ssDNA_overhang_length = ssDNA_overhang_length  # Initial length of ssDNA overhang
          self.ctip_active = ctip_active  # True if CtIP is active
          self.mre11_active = mre11_active  # True if MRE11 is active
          self.exo1_active = exo1_active  # True if EXO1 is active
          self.dna2_active = dna2_active  # True if DNA2 is active

     def start_resection(self):
          # CtIP and MRE11 initiate the resection
          if self.ctip_active and self.mre11_active:
               print("CtIP and MRE11 initiating resection... Trimming 5' ends.")
               self.ssDNA_overhang_length += 10  # CtIP and MRE11 trim some of the 5' ends
          # EXO1 extends the resection
          if self.exo1_active:
               print("EXO1 extending resection... Creating longer ssDNA overhangs.")
               self.ssDNA_overhang_length += 20  # EXO1 extends the ssDNA overhangs
          # DNA2 continues the resection process
          if self.dna2_active:
               print("DNA2 helicase assisting in resection... Further elongating ssDNA overhangs.")
               self.ssDNA_overhang_length += 25  # DNA2 continues to process the overhangs

          print(f"Resection completed. ssDNA overhang length: {self.ssDNA_overhang_length} nucleotides.")

print(model)
print(model.monomers)
print(model.parameters)
print(model.rules)
