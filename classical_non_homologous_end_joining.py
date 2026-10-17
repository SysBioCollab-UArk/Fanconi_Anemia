from pysb import *
from pysb.util import alias_model_components


def create_model_elements(define_observables=True):
    """
    Classical non-homologous end joining (cNHEJ) submodel.

    Design goals
    ------------
    1. Keep the shared DSB definition unchanged:
           Monomer("DSB", ["b"])
       DSB is therefore NOT redefined in this function.

    2. Represent the two physical DNA ends without adding extra sites to DSB.
       Two Ku70/80 molecules are lumped into a single Ku70_80_pair species,
       analogous to the Parp1_dimer strategy used in the MMEJ model.

    3. Match the current cNHEJ pathway diagram:
         A. Ku70/80 end recognition
         B. DNA-PKcs recruitment
         C. Long-range synaptic complex with XRCC4-LIG4 and XLF
         D. Conditional Artemis processing
         E. DNA-PKcs autophosphorylation/rearrangement and release
         F. Short-range synaptic complex
         G. Conditional Pol mu / Pol lambda synthesis
         H. LIG4-mediated ligation
         I. Repaired DNA / complex dissociation

    4. Avoid forcing every DSB through Artemis or polymerase processing.
       The end_state site on Ku70_80_pair is a lumped representation of
       DNA-end chemistry after cNHEJ has engaged the break. It does NOT add
       additional binding sites to the shared DSB monomer.
    """

    # ------------------------------------------------------------------
    # MONOMERS
    # ------------------------------------------------------------------

    # One Ku70/80 heterodimer.
    Monomer("Ku70_80", ["dsb"])

    # Lumped representation of TWO Ku70/80 heterodimers occupying the two
    # physical ends of one DSB. dnapk1 and dnapk2 are bookkeeping labels
    # for the two DNA ends; the numbers do not imply biological asymmetry.
    #
    # stage:
    #   assembling  -> cNHEJ factors are being recruited
    #   long_range  -> long-range synaptic complex
    #   short_range -> ligation-competent short-range complex
    #
    # end_state:
    #   unclassified -> end chemistry has not yet been assigned in this
    #                   lumped model
    #   compatible   -> ends can proceed without Artemis processing
    #   damaged      -> Artemis-dependent processing is required
    #   ready         -> ends are ligation-ready
    #   gap           -> synthesis is still required before ligation
    Monomer("Ku70_80_pair",["dsb", "dnapk1", "dnapk2", "xrcc4", "pol", "stage", "end_state"],
            {"stage": ["assembling", "long_range", "short_range"],
             "end_state": ["unclassified", "compatible", "damaged", "ready", "gap"]})

    # DNA-PK catalytic subunit. "state" lumps the activation/autophosphorylation transition that promotes end
    # accessibility and eventual DNA-PKcs release.
    Monomer("DNA_PKcs", ["ku_pair", "artemis", "state"], {"state": ["inactive", "active"]})

    Monomer("Artemis", ["dnapk"])

    # Alternative cNHEJ polymerases. They compete for the same pol site on the short-range Ku pair and are NOT forced to
    # act sequentially.
    Monomer("PolMu", ["ku_pair"])
    Monomer("PolLambda", ["ku_pair"])

    # Ligation machinery. XRCC4 and LIG4 can preassemble, and XLF joins the scaffold. This allows the ligation machinery
    # to participate in the long-range synaptic complex rather than appearing only at the end.
    Monomer("XRCC4", ["lig4", "xlf", "ku_pair"])
    Monomer("LIG4", ["xrcc4"])
    Monomer("XLF", ["xrcc4"])

    # ------------------------------------------------------------------
    # PARAMETERS / INITIAL CONDITIONS
    # ------------------------------------------------------------------
    # These concentrations and kinetic values are placeholders, consistent with the current project style. They should
    # be fit or revised when experimental constraints are introduced.

    Parameter("Ku70_80_0", 1000)
    Parameter("DNA_PKcs_0", 1000)
    Parameter("Artemis_0", 500)
    Parameter("PolMu_0", 500)
    Parameter("PolLambda_0", 500)
    Parameter("XRCC4_0", 1000)
    Parameter("LIG4_0", 1000)
    Parameter("XLF_0", 1000)

    alias_model_components()

    Initial(Ku70_80(dsb=None), Ku70_80_0)
    Initial(DNA_PKcs(ku_pair=None, artemis=None, state="inactive"), DNA_PKcs_0)
    Initial(Artemis(dnapk=None), Artemis_0)
    Initial(PolMu(ku_pair=None), PolMu_0)
    Initial(PolLambda(ku_pair=None), PolLambda_0)
    Initial(XRCC4(lig4=None, xlf=None, ku_pair=None), XRCC4_0)
    Initial(LIG4(xrcc4=None), LIG4_0)
    Initial(XLF(xrcc4=None), XLF_0)

    # ------------------------------------------------------------------
    # OBSERVABLES
    # ------------------------------------------------------------------
    if define_observables:
        Observable("Ku70_80_monomer", Ku70_80())
        Observable("Ku70_80_pair_total", Ku70_80_pair())

        alias_model_components()

        Expression("Ku70_80_tot", Ku70_80_monomer + 2 * Ku70_80_pair_total)

        Observable("Ku70_80_free", Ku70_80(dsb=None))
        Observable("Ku70_80_single_DSB", DSB(b=1) % Ku70_80(dsb=1))
        Observable("Ku70_80_pair_DSB", DSB(b=1) % Ku70_80_pair(dsb=1))

        Observable("DNA_PKcs_free", DNA_PKcs(ku_pair=None, artemis=None, state="inactive"))
        Observable("DNA_PKcs_bound", DNA_PKcs(ku_pair=ANY))
        Observable("DNA_PKcs_active", DNA_PKcs(state="active"))

        Observable("Artemis_free", Artemis(dnapk=None))
        Observable("Artemis_bound", Artemis(dnapk=ANY))

        Observable("PolMu_free", PolMu(ku_pair=None))
        Observable("PolMu_bound", PolMu(ku_pair=ANY))
        Observable("PolLambda_free", PolLambda(ku_pair=None))
        Observable("PolLambda_bound", PolLambda(ku_pair=ANY))

        Observable("XRCC4_free", XRCC4(lig4=None, xlf=None, ku_pair=None))
        Observable("XRCC4_LIG4", XRCC4(lig4=ANY))
        Observable("XLF_bound", XLF(xrcc4=ANY))

        Observable("cNHEJ_long_range", Ku70_80_pair(stage="long_range"))
        Observable("cNHEJ_short_range", Ku70_80_pair(stage="short_range"))
        Observable("cNHEJ_short_range_gap", Ku70_80_pair(stage="short_range", end_state="gap"))
        Observable("cNHEJ_ligation_ready", Ku70_80_pair(stage="short_range", end_state="ready"))

    # ------------------------------------------------------------------
    # A. DSB RECOGNITION: Ku70/80 BINDS THE DNA ENDS
    # ------------------------------------------------------------------

    Parameter("kf_DSB_Ku70_80", 1)
    Parameter("kr_DSB_Ku70_80", 10)

    alias_model_components()

    # First Ku70/80 binds the one shared PySB DSB binding site.
    Rule("Ku70_80_binds_DSB",
         DSB(b=None) + Ku70_80(dsb=None)
         |
         DSB(b=1) % Ku70_80(dsb=1),
         kf_DSB_Ku70_80, kr_DSB_Ku70_80)

    # A second physical Ku70/80 occupies the other DNA end. Because DSB has
    # only one PySB binding site, the two Ku molecules are replaced by a
    # lumped Ku70_80_pair that remains attached through DSB(b).
    Parameter("kf_second_Ku70_80", 1)
    Parameter("kr_second_Ku70_80", 10)

    alias_model_components()

    Rule("Second_Ku70_80_forms_Ku_pair",
         DSB(b=1) % Ku70_80(dsb=1) + Ku70_80(dsb=None)
         |
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, dnapk2=None, xrcc4=None, pol=None, stage="assembling",
                                 end_state="unclassified"), kf_second_Ku70_80, kr_second_Ku70_80)

    # ------------------------------------------------------------------
    # B. DNA-PK ASSEMBLY: DNA-PKcs BINDS THE TWO Ku-BOUND ENDS
    # ------------------------------------------------------------------

    Parameter("kf_DNA_PKcs_Ku_end1", 1)
    Parameter("kr_DNA_PKcs_Ku_end1", 10)
    Parameter("kf_DNA_PKcs_Ku_end2", 1)
    Parameter("kr_DNA_PKcs_Ku_end2", 10)

    alias_model_components()

    Rule("DNA_PKcs_binds_Ku_pair_end1",
         DNA_PKcs(ku_pair=None, artemis=None, state="inactive") +
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, stage="assembling")
         |
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2,
                                                                                       stage="assembling"),
         kf_DNA_PKcs_Ku_end1, kr_DNA_PKcs_Ku_end1)

    Rule("DNA_PKcs_binds_Ku_pair_end2",
         DNA_PKcs(ku_pair=None, artemis=None, state="inactive") +
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=None,
                                                                                       stage="assembling")
         |
         DNA_PKcs(ku_pair=3, artemis=None, state="inactive") %
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3,
                                                                                       stage="assembling"),
         kf_DNA_PKcs_Ku_end2, kr_DNA_PKcs_Ku_end2)

    # ------------------------------------------------------------------
    # C. LIGATION-SCAFFOLD RECRUITMENT AND LONG-RANGE SYNAPSIS
    # ------------------------------------------------------------------

    # XRCC4 and LIG4 are treated as a preassembling ligation complex.
    Parameter("kf_XRCC4_LIG4", 1)
    Parameter("kr_XRCC4_LIG4", 1)

    alias_model_components()

    Rule("XRCC4_binds_LIG4",
         XRCC4(lig4=None, xlf=None, ku_pair=None) + LIG4(xrcc4=None)
         |
         XRCC4(lig4=1, xlf=None, ku_pair=None) % LIG4(xrcc4=1),
         kf_XRCC4_LIG4, kr_XRCC4_LIG4)

    # The XRCC4-LIG4 complex can be recruited while the Ku/DNA-PK complex is assembling; it is not forced to wait until
    # all end processing is complete.
    Parameter("kf_XRCC4_LIG4_Ku_pair", 1)
    Parameter("kr_XRCC4_LIG4_Ku_pair", 10)

    alias_model_components()

    Rule("XRCC4_LIG4_binds_Ku_pair",
         XRCC4(lig4=4, xlf=None, ku_pair=None) % LIG4(xrcc4=4) +
         DSB(b=1) % Ku70_80_pair(dsb=1, xrcc4=None, stage="assembling")
         |
         DSB(b=1) % Ku70_80_pair(dsb=1, xrcc4=3, stage="assembling") % XRCC4(lig4=4, xlf=None, ku_pair=3) %
         LIG4(xrcc4=4),
         kf_XRCC4_LIG4_Ku_pair, kr_XRCC4_LIG4_Ku_pair)

    # XLF joins the XRCC4-LIG4 scaffold.
    Parameter("kf_XLF_XRCC4", 1)
    Parameter("kr_XLF_XRCC4", 10)

    alias_model_components()

    Rule("XLF_binds_XRCC4_at_DSB",
         XLF(xrcc4=None) + DSB(b=1) % Ku70_80_pair(dsb=1, xrcc4=3, stage="assembling") %
         XRCC4(lig4=4, xlf=None, ku_pair=3) % LIG4(xrcc4=4)
         |
         XLF(xrcc4=5) % DSB(b=1) % Ku70_80_pair(dsb=1, xrcc4=3, stage="assembling") % XRCC4(lig4=4, xlf=5, ku_pair=3) %
         LIG4(xrcc4=4),
         kf_XLF_XRCC4, kr_XLF_XRCC4)

    # Once the two DNA-PKcs molecules and XRCC4-LIG4-XLF scaffold are
    # present, the complex enters the long-range synaptic state.
    Parameter("k_long_range_synapsis", 1)

    alias_model_components()

    Rule("Long_range_synapsis",
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DNA_PKcs(ku_pair=3, artemis=None, state="inactive") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, xrcc4=4, pol=None, stage="assembling",
                                 end_state="unclassified") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) %
         XLF(xrcc4=6)
         >>
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DNA_PKcs(ku_pair=3, artemis=None, state="inactive") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, xrcc4=4, pol=None, stage="long_range",
                                 end_state="unclassified") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) %
         XLF(xrcc4=6),
         k_long_range_synapsis)

    # ------------------------------------------------------------------
    # DNA-PKcs ACTIVATION / AUTOPHOSPHORYLATION
    # ------------------------------------------------------------------

    # This is a lumped complex-level activation step. It avoids modeling
    # individual DNA-PKcs phosphorylation sites while preserving the
    # mechanistically important protected -> processing/release transition.
    Parameter("k_DNA_PKcs_autophosphorylation", 1)

    alias_model_components()

    Rule("DNA_PKcs_autophosphorylation",
         DNA_PKcs(ku_pair=2, artemis=None, state="inactive") % DNA_PKcs(ku_pair=3, artemis=None, state="inactive") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range")
         >>
         DNA_PKcs(ku_pair=2, artemis=None, state="active") % DNA_PKcs(ku_pair=3, artemis=None, state="active") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range"),
         k_DNA_PKcs_autophosphorylation)

    # ------------------------------------------------------------------
    # BRANCH: END STRUCTURE / PROCESSING REQUIREMENT
    # ------------------------------------------------------------------
    # The shared DSB monomer does not encode DNA-end chemistry. Therefore the cNHEJ-specific Ku-pair carries a lumped
    # end_state once long-range synapsis has formed. These two rules represent the fraction of breaks whose ends are
    # already compatible versus those requiring Artemis.
    #
    # Their relative rates can later be fit or replaced by experimentally constrained values.

    Parameter("k_ends_compatible", 1)
    Parameter("k_ends_damaged", 1)

    alias_model_components()

    Rule("Long_range_ends_are_compatible",
         Ku70_80_pair(stage="long_range", end_state="unclassified", dnapk1=ANY, dnapk2=ANY)
         >>
         Ku70_80_pair(stage="long_range", end_state="compatible", dnapk1=ANY, dnapk2=ANY),
         k_ends_compatible)

    Rule("Long_range_ends_require_processing",
         Ku70_80_pair(stage="long_range", end_state="unclassified", dnapk1=ANY, dnapk2=ANY)
         >>
         Ku70_80_pair(stage="long_range", end_state="damaged", dnapk1=ANY, dnapk2=ANY),
         k_ends_damaged)

    # ------------------------------------------------------------------
    # D. CONDITIONAL ARTEMIS PROCESSING
    # ------------------------------------------------------------------

    Parameter("kf_Artemis_DNA_PKcs", 1)
    Parameter("kr_Artemis_DNA_PKcs", 10)

    alias_model_components()

    # One Artemis molecule represents the nuclease-processing event at the
    # synapsed break. It binds an activated DNA-PKcs only when end processing
    # is required.
    Rule("Artemis_binds_active_DNA_PKcs",
         Artemis(dnapk=None) + DNA_PKcs(ku_pair=2, artemis=None, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="damaged")
         |
         Artemis(dnapk=7) % DNA_PKcs(ku_pair=2, artemis=7, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="damaged"),
         kf_Artemis_DNA_PKcs, kr_Artemis_DNA_PKcs)

    # Artemis processing can make the termini directly ligation-ready...
    Parameter("k_Artemis_generates_ready_ends", 1)
    # ...or can leave a short gap that requires polymerase synthesis.
    Parameter("k_Artemis_leaves_gap", 1)

    alias_model_components()

    Rule("Artemis_processing_generates_ready_ends",
         Artemis(dnapk=7) % DNA_PKcs(ku_pair=2, artemis=7, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="damaged")
         >>
         Artemis(dnapk=None) + DNA_PKcs(ku_pair=2, artemis=None, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="ready"),
         k_Artemis_generates_ready_ends)

    Rule("Artemis_processing_leaves_gap",
         Artemis(dnapk=7) % DNA_PKcs(ku_pair=2, artemis=7, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="damaged")
         >>
         Artemis(dnapk=None) + DNA_PKcs(ku_pair=2, artemis=None, state="active") %
         DNA_PKcs(ku_pair=3, artemis=None, state="active") % DSB(b=1) %
         Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, stage="long_range", end_state="gap"),
         k_Artemis_leaves_gap)

    # ------------------------------------------------------------------
    # E/F. DNA-PKcs REARRANGEMENT/RELEASE -> SHORT-RANGE SYNAPSIS
    # ------------------------------------------------------------------
    # These rules lump DNA-PKcs rearrangement and departure into the long-range -> short-range transition. XRCC4-LIG4
    # and XLF remain at the break.

    Parameter("k_DNA_PKcs_release_compatible", 1)
    Parameter("k_DNA_PKcs_release_ready", 1)
    Parameter("k_DNA_PKcs_release_gap", 1)

    alias_model_components()

    Rule("DNA_PKcs_release_from_compatible_ends",
         DNA_PKcs(ku_pair=2, artemis=None, state="active") % DNA_PKcs(ku_pair=3, artemis=None, state="active") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, xrcc4=4, pol=None, stage="long_range",
                                 end_state="compatible") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) %
         XLF(xrcc4=6)
         >>
         DNA_PKcs(ku_pair=None, artemis=None, state="active") + DNA_PKcs(ku_pair=None, artemis=None, state="active") +
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, dnapk2=None, xrcc4=4, pol=None, stage="short_range",
                                 end_state="ready") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6),
         k_DNA_PKcs_release_compatible)

    Rule("DNA_PKcs_release_after_processing_ready",
         DNA_PKcs(ku_pair=2, artemis=None, state="active") % DNA_PKcs(ku_pair=3, artemis=None, state="active") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, xrcc4=4, pol=None, stage="long_range", end_state="ready") %
         XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6)
         >>
         DNA_PKcs(ku_pair=None, artemis=None, state="active") + DNA_PKcs(ku_pair=None, artemis=None, state="active") +
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, dnapk2=None, xrcc4=4, pol=None, stage="short_range",
                                 end_state="ready") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6),
         k_DNA_PKcs_release_ready)

    Rule("DNA_PKcs_release_after_processing_gap",
         DNA_PKcs(ku_pair=2, artemis=None, state="active") % DNA_PKcs(ku_pair=3, artemis=None, state="active") %
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=2, dnapk2=3, xrcc4=4, pol=None, stage="long_range", end_state="gap") %
         XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6)
         >>
         DNA_PKcs(ku_pair=None, artemis=None, state="active") + DNA_PKcs(ku_pair=None, artemis=None, state="active") +
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, dnapk2=None, xrcc4=4, pol=None, stage="short_range",
                                 end_state="gap") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6),
         k_DNA_PKcs_release_gap)

    # Free activated DNA-PKcs is recycled back to the inactive pool.
    Parameter("k_DNA_PKcs_dephosphorylation", 1)

    alias_model_components()

    Rule("DNA_PKcs_dephosphorylation",
         DNA_PKcs(ku_pair=None, artemis=None, state="active")
         >>
         DNA_PKcs(ku_pair=None, artemis=None, state="inactive"),
         k_DNA_PKcs_dephosphorylation)

    # ------------------------------------------------------------------
    # G. CONDITIONAL POLYMERASE STEP
    # ------------------------------------------------------------------
    # Only gap-containing short-range complexes recruit polymerase.
    # Pol mu and Pol lambda compete for the SAME site and are alternatives.

    Parameter("kf_PolMu_short_range", 1)
    Parameter("kr_PolMu_short_range", 10)
    Parameter("kf_PolLambda_short_range", 1)
    Parameter("kr_PolLambda_short_range", 10)
    Parameter("k_PolMu_gap_fill", 1)
    Parameter("k_PolLambda_gap_fill", 1)

    alias_model_components()

    Rule("PolMu_binds_short_range_gap",
         PolMu(ku_pair=None) + DSB(b=1) % Ku70_80_pair(dsb=1, pol=None, stage="short_range", end_state="gap")
         |
         PolMu(ku_pair=7) % DSB(b=1) % Ku70_80_pair(dsb=1, pol=7, stage="short_range", end_state="gap"),
         kf_PolMu_short_range, kr_PolMu_short_range)

    Rule("PolLambda_binds_short_range_gap",
         PolLambda(ku_pair=None) + DSB(b=1) % Ku70_80_pair(dsb=1, pol=None, stage="short_range", end_state="gap")
         |
         PolLambda(ku_pair=7) % DSB(b=1) % Ku70_80_pair(dsb=1, pol=7, stage="short_range", end_state="gap"),
         kf_PolLambda_short_range, kr_PolLambda_short_range)

    Rule("PolMu_fills_gap",
         PolMu(ku_pair=7) % DSB(b=1) % Ku70_80_pair(dsb=1, pol=7, stage="short_range", end_state="gap")
         >>
         PolMu(ku_pair=None) + DSB(b=1) % Ku70_80_pair(dsb=1, pol=None, stage="short_range", end_state="ready"),
         k_PolMu_gap_fill)

    Rule("PolLambda_fills_gap",
         PolLambda(ku_pair=7) % DSB(b=1) % Ku70_80_pair(dsb=1, pol=7, stage="short_range", end_state="gap")
         >>
         PolLambda(ku_pair=None) + DSB(b=1) % Ku70_80_pair(dsb=1, pol=None, stage="short_range", end_state="ready"),
         k_PolLambda_gap_fill)

    # ------------------------------------------------------------------
    # H/I. LIG4-MEDIATED LIGATION AND REPAIR-COMPLEX DISSOCIATION
    # ------------------------------------------------------------------

    Parameter("k_LIG4_ligates_DSB", 1)

    alias_model_components()

    # Successful ligation consumes the DSB species. The lumped Ku pair is converted back into two free Ku70/80
    # heterodimers, while XRCC4, LIG4, and XLF are recycled.
    Rule("LIG4_ligates_DSB",
         DSB(b=1) % Ku70_80_pair(dsb=1, dnapk1=None, dnapk2=None, xrcc4=4, pol=None, stage="short_range",
                                 end_state="ready") % XRCC4(lig4=5, xlf=6, ku_pair=4) % LIG4(xrcc4=5) % XLF(xrcc4=6)
         >>
         Ku70_80(dsb=None) + Ku70_80(dsb=None) + XRCC4(lig4=None, xlf=None, ku_pair=None) + LIG4(xrcc4=None) +
         XLF(xrcc4=None),
         k_LIG4_ligates_DSB)


if __name__ == "__main__":
    # Standalone test harness. Shared components are defined here only for testing this module by itself. In
    # FA_merged_model.py, DSB should remain defined centrally as Monomer("DSB", ["b"]).

    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()

    Monomer("DSB", ["b"])
    Parameter("DSB_0", 100)
    Initial(DSB(b=None), DSB_0)
    Observable("DSB_tot", DSB())
    Observable("DSB_free", DSB(b=None))

    create_model_elements(define_observables=True)

    tspan = np.linspace(0, 20, 2001)
    sim = ScipyOdeSimulator(model, tspan=tspan, verbose=True)
    output = sim.run()

    print(model.observables)

    # DSB dynamics
    plt.figure(constrained_layout=True)
    for obs in [DSB_tot, DSB_free]:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.legend(loc="best")

    # cNHEJ progression
    plt.figure(constrained_layout=True)
    for obs in [Ku70_80_single_DSB, Ku70_80_pair_DSB, cNHEJ_long_range, cNHEJ_short_range, cNHEJ_short_range_gap,
                cNHEJ_ligation_ready]:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.legend(loc="best")

    # Processing proteins
    plt.figure(constrained_layout=True)
    for obs in [DNA_PKcs_bound, DNA_PKcs_active, Artemis_bound, PolMu_bound, PolLambda_bound]:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.legend(loc="best")

    plt.show()
