module Instrument

  export rosetta_instruments

immutable instrument
  nPixelsX::Int
  nPixelsY::Int
  phiX::Float64
  phiY::Float64
  frame::AbstractString
end

osiris_wac = instrument(1024, 1024, 6.0, 6.0, "ROS_OSIRIS_WAC")
osiris_nac = instrument(1024, 1024, 6.0, 6.0, "ROS_OSIRIS_NAC")
alice = instrument(19, 1, 5.852/2.0, 0.0, "ROS_ALICE")
miro = instrument(1, 1, 0.0, 0.0, "ROS_MIRO_MM")
virtis_m = instrument(256, 256, 3.6669/2.0, 3.6669/2.0, "ROS_VIRTIS-M")
debug = instrument(100,100,1.5,1.5, "ROS_NADIR")

rosetta_instruments = Dict()
rosetta_instruments["OSIRIS_WAC"] = osiris_wac
rosetta_instruments["OSIRIS_NAC"] = osiris_nac
rosetta_instruments["ALICE"] = alice
rosetta_instruments["MIRO"] = miro
rosetta_instruments["VIRTIS_M"] = virtis_m
rosetta_instruments["DEBUG"] = debug

end
