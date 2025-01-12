import numpy as np
import math

def score2(fuel_weight, time):
    Max_feul_weith_time_ratio = 3 / 90 # 얼마로 나오려나
    return 1 + ((fuel_weight / time) / Max_feul_weith_time_ratio)

def score3(laps_flown, X1_weight):
    Max_laps_flown_X1_weight_ratio = 5 + 2.5 / 0.3 # 얼마로 나오려나
    bonus_box_score = 2.5
    return 2 + ((laps_flown + bonus_box_score / X1_weight) / Max_laps_flown_X1_weight_ratio)