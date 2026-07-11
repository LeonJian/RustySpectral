# Supported Platforms

RustySpectral supports the same satellite platforms and sensors as pyspectral
through an identical INSTRUMENTS mapping (56 entries).

## Platform → Sensor Mapping

### Geostationary

| Platform | Sensor | RSR Filename |
|----------|--------|-------------|
| Meteosat-8 | seviri | `rsr_seviri_Meteosat-8.h5` |
| Meteosat-9 | seviri | `rsr_seviri_Meteosat-9.h5` |
| Meteosat-10 | seviri | `rsr_seviri_Meteosat-10.h5` |
| Meteosat-11 | seviri | `rsr_seviri_Meteosat-11.h5` |
| GOES-16 | abi | `rsr_abi_GOES-16.h5` |
| GOES-17 | abi | `rsr_abi_GOES-17.h5` |
| GOES-18 | abi | `rsr_abi_GOES-18.h5` |
| GOES-19 | abi | `rsr_abi_GOES-19.h5` |
| Himawari-8 | ahi | `rsr_ahi_Himawari-8.h5` |
| Himawari-9 | ahi | `rsr_ahi_Himawari-9.h5` |
| GEO-KOMPSAT-2A | ami | `rsr_ami_GEO-KOMPSAT-2A.h5` |
| GEO-KOMPSAT-2B | goci-2 | `rsr_goci-2_GEO-KOMPSAT-2B.h5` |
| FY-4A | agri | `rsr_agri_FY-4A.h5` |
| FY-4B | agri / ghi | `rsr_agri_FY-4B.h5` |
| Meteosat-12 | fci | `rsr_fci_Meteosat-12.h5` |
| MTG-I1 | fci | `rsr_fci_MTG-I1.h5` |
| Electro-L-N2 | msu-gs | `rsr_msu-gs_Electro-L-N2.h5` |
| Arctica-M-N1 | msu-gsa | `rsr_msu-gsa_Arctica-M-N1.h5` |

### Polar-Orbiting (LEO)

| Platform | Sensor | RSR Filename |
|----------|--------|-------------|
| TIROS-N | avhrr/1 | `rsr_avhrr1_TIROS-N.h5` |
| NOAA-6 | avhrr/1 | `rsr_avhrr1_NOAA-6.h5` |
| NOAA-8 | avhrr/1 | `rsr_avhrr1_NOAA-8.h5` |
| NOAA-10 | avhrr/1 | `rsr_avhrr1_NOAA-10.h5` |
| NOAA-7 | avhrr/2 | `rsr_avhrr2_NOAA-7.h5` |
| NOAA-9 | avhrr/2 | `rsr_avhrr2_NOAA-9.h5` |
| NOAA-11 | avhrr/2 | `rsr_avhrr2_NOAA-11.h5` |
| NOAA-12 | avhrr/2 | `rsr_avhrr2_NOAA-12.h5` |
| NOAA-14 | avhrr/2 | `rsr_avhrr2_NOAA-14.h5` |
| NOAA-15 | avhrr/3 | `rsr_avhrr3_NOAA-15.h5` |
| NOAA-16 | avhrr/3 | `rsr_avhrr3_NOAA-16.h5` |
| NOAA-17 | avhrr/3 | `rsr_avhrr3_NOAA-17.h5` |
| NOAA-18 | avhrr/3 | `rsr_avhrr3_NOAA-18.h5` |
| NOAA-19 | avhrr/3 | `rsr_avhrr3_NOAA-19.h5` |
| Metop-A | avhrr/3 | `rsr_avhrr3_Metop-A.h5` |
| Metop-B | avhrr/3 | `rsr_avhrr3_Metop-B.h5` |
| Metop-C | avhrr/3 | `rsr_avhrr3_Metop-C.h5` |
| EOS-Aqua | modis | `rsr_modis_EOS-Aqua.h5` |
| EOS-Terra | modis | `rsr_modis_EOS-Terra.h5` |
| Aqua | modis | `rsr_modis_EOS-Aqua.h5` |
| Terra | modis | `rsr_modis_EOS-Terra.h5` |
| Suomi-NPP | viirs | `rsr_viirs_Suomi-NPP.h5` |
| NOAA-20 | viirs | `rsr_viirs_NOAA-20.h5` |
| NOAA-21 | viirs | `rsr_viirs_NOAA-21.h5` |
| Metop-SG-A1 | metimage | `rsr_metimage_Metop-SG-A1.h5` |
| FY-3A | virr / mersi-1 | `rsr_virr_FY-3A.h5` |
| FY-3B | virr / mersi-1 | `rsr_virr_FY-3B.h5` |
| FY-3C | virr / mersi-1 | `rsr_virr_FY-3C.h5` |
| FY-3D | mersi-2 | `rsr_mersi-2_FY-3D.h5` |
| FY-3F | mersi-3 | `rsr_mersi-3_FY-3F.h5` |
| FY-3G | mersi-rm | `rsr_mersi-rm_FY-3G.h5` |

### Ocean Color

| Platform | Sensor | RSR Filename |
|----------|--------|-------------|
| Sentinel-3A | olci / slstr | `rsr_olci_Sentinel-3A.h5` |
| Sentinel-3B | olci / slstr | `rsr_olci_Sentinel-3B.h5` |
| HY-1C | cocts | `rsr_cocts_HY-1C.h5` |
| Envisat | aatsr | `rsr_aatsr_Envisat.h5` |

### Land Imaging

| Platform | Sensor | RSR Filename |
|----------|--------|-------------|
| Sentinel-2A | msi | `rsr_msi_Sentinel-2A.h5` |
| Sentinel-2B | msi | `rsr_msi_Sentinel-2B.h5` |
| Sentinel-2C | msi | `rsr_msi_Sentinel-2C.h5` |
| Landsat-8 | oli_tirs | `rsr_oli_tirs_Landsat-8.h5` |
| Landsat-9 | oli_tirs | `rsr_oli_tirs_Landsat-9.h5` |

### Deep Space

| Platform | Sensor | RSR Filename |
|----------|--------|-------------|
| DSCOVR | epic | `rsr_epic_DSCOVR.h5` |

## Band Name Dictionaries

Band name lookups are available for the following 12 instrument types:

| Key | Instrument | Band Count |
|-----|-----------|------------|
| `generic` | Generic (all instruments) | 72 |
| `modis` | Terra/Aqua MODIS | 36 |
| `seviri` | Meteosat SEVIRI | 12 |
| `viirs` | Suomi-NPP/NOAA-20 VIIRS | 14 |
| `avhrr3` | AVHRR/3 | 6 |
| `abi` | GOES-R ABI | 16 |
| `agri` | FY-4 AGRI | 14 |
| `ahi` | Himawari AHI | 16 |
| `ami` | GEO-KOMPSAT-2A AMI | 16 |
| `fci` | MTG FCI | 16 |
| `slstr` | Sentinel-3 SLSTR | 11 |
| `vii` | EPS-SG VII | 20 |

### SEVIRI Wave Regression Parameters

The following regression coefficients are built in for SEVIRI radiance ↔
brightness temperature conversion:

| Band | Meteosat-8 (ν, α, β) | Meteosat-9 (ν, α, β) |
|------|----------------------|----------------------|
| IR3.9 | 2567.330, 0.9956, 3.410 | 2568.832, 0.9954, 3.438 |
| WV6.2 | 1598.103, 0.9962, 2.218 | 1600.548, 0.9963, 2.185 |
| WV7.3 | 1362.081, 0.9991, 0.478 | 1360.330, 0.9991, 0.470 |
| IR8.7 | 1149.069, 0.9996, 0.179 | 1148.620, 0.9996, 0.179 |
| IR9.7 | 1034.343, 0.9999, 0.060 | 1035.289, 0.9999, 0.056 |
| IR10.8 | 930.647, 0.9983, 0.625 | 931.700, 0.9983, 0.640 |
| IR12.0 | 839.660, 0.9988, 0.397 | 836.445, 0.9988, 0.408 |
| IR13.4 | 752.387, 0.9981, 0.578 | 751.792, 0.9981, 0.561 |

ν in cm⁻¹, α dimensionless, β in K.

## Atmosphere Types

| Identifier (underscore) | Identifier (Python-style) |
|------------------------|--------------------------|
| `subarctic_summer` | `subarctic summer` |
| `subarctic_winter` | `subarctic winter` |
| `midlatitude_summer` | `midlatitude summer` |
| `midlatitude_winter` | `midlatitude winter` |
| `tropical` | `tropical` |
| `us_standard` | `us-standard` |

## Aerosol Types

`antarctic_aerosol`, `continental_average_aerosol`, `continental_clean_aerosol`,
`continental_polluted_aerosol`, `desert_aerosol`, `marine_clean_aerosol`,
`marine_polluted_aerosol`, `marine_tropical_aerosol`, `rayleigh_only`,
`rural_aerosol`, `urban_aerosol`
