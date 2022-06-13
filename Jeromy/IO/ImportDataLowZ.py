import Jeromy
import Jeromy.IO.FileImporter as FI

# B8_CC      = FI.ImportFile("Boron8 CC", Jeromy.__jeromy_data__, "B8_CC_smooth_lowZ.txt" )
# B8_ES      = FI.ImportFile("Boron8 ES", Jeromy.__jeromy_data__, "B8_ES_smooth_lowZ.txt" )
# B8         = FI.ImportFile("Boron8"   , Jeromy.__jeromy_data__, "B8_ES_smooth_lowZ.txt" )
# B8_map     = FI.ImportFile.combine_inputs([B8_CC, B8_ES])
# B8.map     = B8_map
# B8.yValues = [B8.map[i] for i in B8.map]


# HEP_CC      = FI.ImportFile("HEP CC"   , Jeromy.__jeromy_data__, "HEP_CC_smooth_lowZ.txt")
# HEP_ES      = FI.ImportFile("HEP ES"   , Jeromy.__jeromy_data__, "HEP_ES_smooth_lowZ.txt")
# HEP         = FI.ImportFile("HEP"      , Jeromy.__jeromy_data__, "B8_ES_smooth_lowZ.txt" )
# HEP_map     = FI.ImportFile.combine_inputs([HEP_CC, HEP_ES])
# HEP.map     = HEP_map
# HEP.yValues = [HEP.map[i] for i in HEP.map]

# O15         = FI.ImportFile("O15 "     , Jeromy.__jeromy_data__, "O15_smooth_lowZ.txt"   )
# F17         = FI.ImportFile("F17 "     , Jeromy.__jeromy_data__, "F17_smooth_lowZ.txt"   )
# CNO         = FI.ImportFile("CNO "     , Jeromy.__jeromy_data__, "F17_smooth_lowZ.txt"   )
# CNO_map     = FI.ImportFile.combine_inputs([O15, F17])
# CNO.map     = CNO_map
# CNO.yValues = [CNO.map[i] for i in CNO.map]


# Radon                 = FI.ImportFile("Radon"        , Jeromy.__jeromy_data__, "radon_100kt.txt"  )
# Neutron               = FI.ImportFile("Neutron"      , Jeromy.__jeromy_data__, "neutron_100kt.txt")
# Ar42                  = FI.ImportFile("Ar42"         , Jeromy.__jeromy_data__, "ar42_100kt.txt"   )
# Radiologicals         = FI.ImportFile("Radiologicals", Jeromy.__jeromy_data__, "radon_100kt.txt"  )
# Radiologicals_map     = FI.ImportFile.combine_inputs([Radon, Neutron, Ar42])
# Radiologicals.map     = Radiologicals_map
# Radiologicals.yValues = [Radiologicals.map[i] for i in Radiologicals.map]

# Backgrounds         = FI.ImportFile("Backgrounds", Jeromy.__jeromy_data__, "radon_100kt.txt"  )
# Backgrounds_map     = FI.ImportFile.combine_inputs([Radon, Neutron, Ar42, B8, HEP])
# Backgrounds.map     = Backgrounds_map
# Backgrounds.yValues = [Backgrounds.map[i] for i in Backgrounds.map]

# TotalSignal         = FI.ImportFile("TotalSignal", Jeromy.__jeromy_data__, "radon_100kt.txt"  )
# TotalSignal_map     = FI.ImportFile.combine_inputs([Radon, Neutron, Ar42, B8, HEP, O15, F17])
# TotalSignal.map     = TotalSignal_map
# TotalSignal.yValues = [TotalSignal.map[i] for i in TotalSignal.map]

B8_CC      = FI.ImportFile("Boron8 CC", Jeromy.__jeromy_data__, "B8_CC_smooth_lowZ.txt" )
B8_ES      = FI.ImportFile("Boron8 ES", Jeromy.__jeromy_data__, "B8_ES_smooth_lowZ.txt" )
B8         = FI.ImportFile(name="Boron8")
B8.combine_inputs([B8_CC, B8_ES])
# B8.map     = B8_map
# B8.yValues = [B8.map[i] for i in B8.map]


HEP_CC      = FI.ImportFile("HEP CC"   , Jeromy.__jeromy_data__, "HEP_CC_smooth_lowZ.txt")
HEP_ES      = FI.ImportFile("HEP ES"   , Jeromy.__jeromy_data__, "HEP_ES_smooth_lowZ.txt")
HEP         = FI.ImportFile("HEP")
HEP.combine_inputs([HEP_CC, HEP_ES])
# HEP.map     = HEP_map
# HEP.yValues = [HEP.map[i] for i in HEP.map]

O15         = FI.ImportFile("O15 "     , Jeromy.__jeromy_data__, "O15_smooth_lowZ.txt"   )
F17         = FI.ImportFile("F17 "     , Jeromy.__jeromy_data__, "F17_smooth_lowZ.txt"   )
CNO         = FI.ImportFile("CNO ")
CNO.combine_inputs([O15, F17])



Radon                 = FI.ImportFile("Radon"        , Jeromy.__jeromy_data__, "radon_100kt.txt"  )
Neutron               = FI.ImportFile("Neutron"      , Jeromy.__jeromy_data__, "neutron_100kt.txt")
Ar42                  = FI.ImportFile("Ar42"         , Jeromy.__jeromy_data__, "ar42_100kt.txt"   )
Radiologicals         = FI.ImportFile("Radiologicals")
Radiologicals.combine_inputs([Radon, Neutron, Ar42])


Backgrounds         = FI.ImportFile("Backgrounds")
Backgrounds.combine_inputs([Radon, Neutron, Ar42, B8, HEP])

TotalSignal         = FI.ImportFile("TotalSignal")
TotalSignal.combine_inputs([Radon, Neutron, Ar42, B8, HEP, O15, F17])
# TotalSignal.map     = TotalSignal_map
# TotalSignal.yValues = [TotalSignal.map[i] for i in TotalSignal.map]

