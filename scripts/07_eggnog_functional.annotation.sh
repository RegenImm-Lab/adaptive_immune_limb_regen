#eggnog used for functional annotation, version info below
emapper.py --version
#emapper-2.1.13 / Expected eggNOG DB version: 5.0.2 / Installed eggNOG DB version: 5.0.2 / Diamond version found: diamond version 2.0.15 / MMseqs2 version found: 16.747c6 / Compatible novel families DB version: 1.0.1

#perform functional annotation of proteome, keep in mind this is annotating each isoform so one needs to still collapse this using the provided 20251110_collapse_eggnog.R
emapper.py --cpu 48 -m mmseqs --itype proteins -i ${WD}GCF_040938575.1_UKY_AmexF1_1_protein.faa  -o ${OUT}axolotl_eggnog
