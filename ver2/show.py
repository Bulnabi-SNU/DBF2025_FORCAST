import argparse
from vsp_analysis import loadAnalysisResults, visualize_results
import pandas as pd
from mission_analysis import MissionAnalyzer, visualize_mission
from models import MissionParameters
from config import PresetValues, PropulsionSpecs
import numpy as np

def get_result_by_id(resultID:str, csvPath: str="data/total_results.csv")->pd.DataFrame:
    resultID_df = pd.read_csv(csvPath, sep='|',encoding='utf-8')
    resultID_df = resultID_df[resultID_df['resultID'] == resultID]
    return resultID_df
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="module that displays the result screen which user wants.")
    subparsers = parser.add_subparsers(dest="main_command", required=True)
    show_parser = subparsers.add_parser("show", help="Show results for specific resultID.")
    
    show_subparsers = show_parser.add_subparsers(dest="type", required=True)
    
    show_aircraft_parser = show_subparsers.add_parser("aircraft", help="Show aircraft analysis results.")
    show_aircraft_parser.add_argument("hashVal", type=int, help="Enter the aircraft hash which you want to check.")
    
    show_mission_parser = show_subparsers.add_parser("mission2", help="Show aircraft analysis results.")
    show_mission_parser.add_argument("resultID", type=str, help="Enter the resultID which you want to check.")
    
    show_mission_parser = show_subparsers.add_parser("mission3", help="Show aircraft analysis results.")
    show_mission_parser.add_argument("resultID", type=str, help="Enter the resultID which you want to check.")
    
    
    args = parser.parse_args()
    if args.main_command == "show":
        if args.type == "aircraft":
            aircraft_result = loadAnalysisResults(args.hashVal)
            visualize_results(aircraft_result)
            
        elif args.type == "mission2":
            resultID = "'" + args.resultID + "'"
            resultID_df = get_result_by_id(resultID)
            hashVal = resultID_df['hash']  
            aircraft = loadAnalysisResults(hashVal.iloc[0])     
            param2 = MissionParameters(
                            max_battery_capacity = resultID_df['max_battery_capacity'].iloc[0],
                            throttle_takeoff = 0.9,              # Fixed
                            throttle_climb = resultID_df['mission2_throttle_climb'].iloc[0],
                            throttle_level = resultID_df['mission2_throttle_level'].iloc[0],
                            throttle_turn = resultID_df['mission2_throttle_turn'].iloc[0],                
                            max_climb_angle = 40,                # Fixed
                            max_speed= 40,                       # Fixed
                            max_load_factor = 4.0,               # Fixed
                            h_flap_transition = 5                # Fixed
            )
            
            presetValues = PresetValues(
                            m_x1= resultID_df['m_x1'].iloc[0],
                            x1_flight_time= resultID_df['x1_flight_time'].iloc[0],

                            max_battery_capacity= resultID_df['max_battery_capacity'].iloc[0],
                            min_battery_voltage= resultID_df['min_battery_voltage'].iloc[0],

                            Thrust_max= resultID_df['Thrust_max'].iloc[0],
                            propulsion_efficiency= resultID_df['propulsion_efficiency'].iloc[0],
                            score_weight_ratio= resultID_df['score_weight_ratio'].iloc[0]               
                
            )
              
            propulsionSpecs = PropulsionSpecs(
            propeller_data_path = "data/Propeller10x6E.csv",
            battery_data_path = "data/Maxamps_2250mAh_6S.csv",
            Kv = 109.91,
            R = 0.062,
            max_current = 60,
            max_power = 1332    
            )  
              
                
            missionAnalyzer2 = MissionAnalyzer(aircraft,param2,presetValues, propulsionSpecs)
            missionAnalyzer2.run_mission2()
            visualize_mission(missionAnalyzer2.stateLog)  

        elif args.type == "mission3":
            resultID = "'" + args.resultID + "'"
            resultID_df = get_result_by_id(resultID)
            hashVal = resultID_df['hash']  
            aircraft = loadAnalysisResults(hashVal.iloc[0])     

            param3 = MissionParameters(
                            max_battery_capacity = resultID_df['max_battery_capacity'].iloc[0],
                            throttle_takeoff = 0.9,              # Fixed
                            throttle_climb = resultID_df['mission3_throttle_climb'].iloc[0],
                            throttle_level = resultID_df['mission3_throttle_level'].iloc[0],
                            throttle_turn = resultID_df['mission3_throttle_turn'].iloc[0],                
                            max_climb_angle = 40,                # Fixed
                            max_speed= 40,                       # Fixed
                            max_load_factor = 4.0,               # Fixed
                            h_flap_transition = 5                # Fixed
            )
            
            presetValues = PresetValues(
                            m_x1= resultID_df['m_x1'].iloc[0],
                            x1_flight_time= resultID_df['x1_flight_time'].iloc[0],

                            max_battery_capacity= resultID_df['max_battery_capacity'].iloc[0],
                            min_battery_voltage= resultID_df['min_battery_voltage'].iloc[0],

                            Thrust_max= resultID_df['Thrust_max'].iloc[0],
                            propulsion_efficiency= resultID_df['propulsion_efficiency'].iloc[0],
                            score_weight_ratio= resultID_df['score_weight_ratio'].iloc[0]               
                
            )
             
            propulsionSpecs = PropulsionSpecs(
            propeller_data_path = "data/Propeller10x6E.csv",
            battery_data_path = "data/Maxamps_2250mAh_6S.csv",
            Kv = 109.91,
            R = 0.062,
            max_current = 60,
            max_power = 1332    
            )  
                     
            missionAnalyzer3 = MissionAnalyzer(aircraft,param3,presetValues,propulsionSpecs)
            missionAnalyzer3.run_mission3()
            visualize_mission(missionAnalyzer3.stateLog) 
