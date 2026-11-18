# Third-party components and data

| Component | Use | Licence | Distribution |
|---|---|---|---|
| PyTheranostics 0.1 (Uribe C, QURIT, BC Cancer Research Institute; https://github.com/qurit/PyTheranostics) | 177Lu soft-tissue voxel S-value kernel `doodle/kernels/DK_softtissueICRP_mGyperdecay.img`, resampled to 4.42 mm by step p04 (Supplementary Methods S3) | Apache License 2.0 | Redistributed unmodified, with the licence, in the derived-data record (`derived/kernels/`), because version 0.1 is no longer listed on PyPI and later versions do not contain the kernel. SHA-256 533be3195abccfcdd4ebe8fe0bb8e969539c537f52df2b3e62233d89dbba1022 |
| MRTDosimetry dataset (Robinson AP, Bardiès M, Denis-Bacelar AM, Lassmann M, Tipping J, Tran-Gia J; OSF, doi:10.17605/OSF.IO/69NGE) | measured SPECT/CT, NEMA images, ground-truth concentrations | CC BY 4.0 | Not redistributed (download from OSF; see data/README.md) |
| OpenDose3D (3D Slicer extension) | absorbed-dose rates of the GT maps and of the measured SPECT (outputs in the derived-data record) | see the OpenDose3D distribution | Not included |
