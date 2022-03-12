probe_collection = []

# from compare_DENQ_probe import CompareDENQDensities
# probe_collection.append(CompareDENQDensities)

from ion_screen_probe import AmIAnIonML
# not used: AmIAnIon (original), runs VERY SLOWLY on larger models
probe_collection.append(AmIAnIonML)

from ptxm_probe import IsModifiedNucleotide
probe_collection.append(IsModifiedNucleotide)

# from map_CC_probe import CryoMapModelCC, XrayMapModelCC, NeutronMapModelCC
# probe_collection.extend([CryoMapModelCC, XrayMapModelCC, NeutronMapModelCC])