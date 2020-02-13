# hfe_cleaner.py
# simple script to apply corrections to apollo hfe data

from hfe_corrections import *
from hfe_utils import *

print("Ingesting NSSDC data.")
data = load_hfe(pathname=".")

print("Ingesting 2018 Nagihara data.")
nagihara_data = ingest_nagihara_2018(
    spreadsheet_path="./nagihara/nagihara_jgr_planets.xlsx"
)

print("Ingesting 2019 Nagihara data.")
nagihara_data = ingest_nagihara_2019(nagihara_data, pathname="./nagihara/pds")

data = {**data, **nagihara_data}

print("Flagging missing data.")
flag_missing_hfe(data)

print("Applying per-file flags and corrections.")
a15p1_1_cleanup(data)
a15_1975p1_1_cleanup(data)
a15p1_2_cleanup(data)
a15_1975_p1_2_cleanup(data)
a15p1_3_cleanup(data)
a15p1_4_cleanup(data)
a15p1_5_cleanup(data)
a15p2_1_cleanup(data)
a15p2_2_cleanup(data)
a15_1975p2_2_cleanup(data)
a15p2_3_cleanup(data)
a15p2_4_cleanup(data)
a15p2_5_cleanup(data)
a17p1_1_cleanup(data)
a17_1975p1_1_cleanup(data)
a17_1976p1_1_cleanup(data)
a17_1977p1_1_cleanup(data)
a17p1_2_cleanup(data)
a17p1_3_cleanup(data)
a17p1_4_cleanup(data)
a17p1_5_cleanup(data)
a17p2_1_cleanup(data)
a17_1975p2_1_cleanup(data)
a17_1976p2_1_cleanup(data)
a17_1977p2_1_cleanup(data)
a17p2_2_cleanup(data)
a17p2_3_cleanup(data)
a17p2_4_cleanup(data)
a17p2_5_cleanup(data)

print("Sorting by time.")
manage_disordered_hfe(data)

print("Writing set 1 ('clean').")

os.mkdir("./data_local")

write_clean_hfe(data, "./data_local/clean")

print("Standardizing time units and discarding undesirable ranges.")
polish_hfe(data)

print("Computing explicit temperature values at each thermometer.")
data_split = split_to_thermometers(data)

print("Concatenating Nagihara and NSSDC data.")
join_hfe_years(data_split)

print("Writing set 2 ('split').")
write_split_hfe(data_split, "./data_local/split")

print("Assigning depth values and concatenating all files per probe.")
data_depth = combine_with_depth(data_split)

print("Writing set 3 ('depth').")
write_deep_hfe(data_depth, "./data_local/depth")

print("Done.")
