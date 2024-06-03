
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20210226_20210325_v001.bsp --output geminids/data/SPICE/spp_recon_20210226_20210325_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20210101_20210226_v001.bsp --output geminids/data/SPICE/spp_recon_20210101_20210226_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20200802_20201016_v001.bsp --output geminids/data/SPICE/spp_recon_20200802_20201016_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20201016_20210101_v001.bsp --output geminids/data/SPICE/spp_recon_20201016_20210101_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20200705_20200802_v001.bsp --output geminids/data/SPICE/spp_recon_20200705_20200802_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20200505_20200705_v001.bsp --output geminids/data/SPICE/spp_recon_20200505_20200705_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20200301_20200505_v001.bsp --output geminids/data/SPICE/spp_recon_20200301_20200505_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20200101_20200301_v001.bsp --output geminids/data/SPICE/spp_recon_20200101_20200301_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20190914_20200101_v001.bsp --output geminids/data/SPICE/spp_recon_20190914_20200101_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20190416_20190914_v001.bsp --output geminids/data/SPICE/spp_recon_20190416_20190914_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20190120_20190416_v001.bsp --output geminids/data/SPICE/spp_recon_20190120_20190416_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20181008_20190120_v001.bsp --output geminids/data/SPICE/spp_recon_20181008_20190120_v001.bsp
# curl https://sppgway.jhuapl.edu/MOC/reconstructed_ephemeris/2018/spp_recon_20180812_20181008_v001.bsp --output geminids/data/SPICE/spp_recon_20180812_20181008_v001.bsp
curl https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls --output geminids/data/SPICE/naif0012.tls
curl https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_part-1.bsp --output geminids/data/SPICE/de431_part-1.bsp
curl https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_part-2.bsp --output geminids/data/SPICE/de431_part-2.bsp
cd geminids/data/SPICE
wget https://www.dropbox.com/scl/fi/ayk4384rjqr9loa7g6vcu/small_kernels.zip\?rlkey=iei7x40aelfo5cuv2cqm0xygx\&dl=1
unzip small_kernels.zip\?rlkey=iei7x40aelfo5cuv2cqm0xygx\&dl=1
rm small_kernels.zip\?rlkey=iei7x40aelfo5cuv2cqm0xygx\&dl=1
cd ../../../data
wget https://www.dropbox.com/scl/fi/orzkj57f4n16hxraalmld/orig_orbit.npy\?rlkey=yxfyhxyvu7cn55e27g85oebgz\&dl=1
mv orig_orbit.npy\?rlkey=yxfyhxyvu7cn55e27g85oebgz\&dl=1 orig_orbit.npy