KPL/MK

   This is the meta-kernel used for the 2021 Geminids Model 
   by Wolf Cukier

   The names and contents of the kernels referenced by this
   meta-kernel are as follows:

   File name                   Contents
   --------------------------  -----------------------------
   naif0008.tls                Generic LSK
   de431.bsp   (part 1 and 2)  Lunar/Planetary Ephemeris SCK
   3003200.bsp                 Phaethon 3200 SCK (via Horizons)
   psp_pedict.bsp              PSP data
   spp_recon*.bsp              More specific pop data
   destiny.bsp                 DESTINY+ Trajectory


   \begindata
   PATH_VALUES = ( 'geminids/data/SPICE' )
   
   PATH_SYMBOLS = ( 'PTH')
   
   KERNELS_TO_LOAD = ( '$PTH/naif0008.tls',
                       '$PTH/de431_part-2.bsp',
                       '$PTH/de431_part-1.bsp',
                       '$PTH/2003200.bsp',
                       '$PTH/psp_predict.bsp',
                       '$PTH/destiny.bsp',
                       '$PTH/spp_recon_20210226_20210325_v001.bsp',
                       '$PTH/spp_recon_20210101_20210226_v001.bsp',
                       '$PTH/spp_recon_20200802_20201016_v001.bsp',
                       '$PTH/spp_recon_20201016_20210101_v001.bsp',
                       '$PTH/spp_recon_20200705_20200802_v001.bsp',
                       '$PTH/spp_recon_20200505_20200705_v001.bsp',
                       '$PTH/spp_recon_20200301_20200505_v001.bsp',
                       '$PTH/spp_recon_20200101_20200301_v001.bsp',
                       '$PTH/spp_recon_20190914_20200101_v001.bsp',
                       '$PTH/spp_recon_20190416_20190914_v001.bsp',
                       '$PTH/spp_recon_20190120_20190416_v001.bsp',
                       '$PTH/spp_recon_20181008_20190120_v001.bsp',
                       '$PTH/spp_recon_20180812_20181008_v001.bsp' )
   \begintext