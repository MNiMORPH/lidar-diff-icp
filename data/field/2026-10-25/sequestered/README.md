# Sequestered: points corrected from the false base B1

`rtk_points_B1_false_base.csv` holds the 15 points (13:20-13:30 CDT, 25 October 2026) whose RTK corrections came
from a student's receiver running as a base, written as exported by `analysis/field_survey/compile_rtk_2026-10-25.py`
with B1's self-reported position alongside. B1's position was its own autonomous average and is wrong by meters
(its frame sits 2.58 m below the 2021 lidar; see the README one level up, "Registration or vegetation?"), so these
rows carry no tie to the DR9925 monument or to the CORS network and must not be joined to `rtk_points_all.csv`.
Held in case a use appears (Andy, 2026-11-20); nothing in the analysis reads this folder.
