probe_collection = []

from ion_screen_probe import AmIAnIonML
probe_collection.append(AmIAnIonML)

from ptxm_probe import IsModifiedNucleotide
probe_collection.append(IsModifiedNucleotide)

# from map_CC_probe import CryoMapModelCC, XrayMapModelCC, NeutronMapModelCC
# probe_collection.extend([CryoMapModelCC, XrayMapModelCC, NeutronMapModelCC])