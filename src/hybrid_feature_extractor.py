import math
import os
import os.path
from collections import defaultdict
from copy import deepcopy
from typing import Union

import numpy as np
import pandas as pd
import scipy

from simba.utils.read_write import read_df, get_fn_ext, write_df, find_core_cnt
from simba.mixins.config_reader import ConfigReader
from simba.mixins.feature_extraction_mixin import FeatureExtractionMixin
from simba.mixins.geometry_mixin import GeometryMixin
from simba.mixins.abstract_classes import AbstractFeatureExtraction
from simba.utils.errors import NoFilesFoundError
from simba.utils.printing import SimbaTimer, stdout_success
import argparse

class ExtractFeaturesFrom5bps(ConfigReader,
                              FeatureExtractionMixin,
                              GeometryMixin,
                              AbstractFeatureExtraction
                              ):
    
    def __init__(self,
                 config_path: str):
        
        FeatureExtractionMixin.__init__(self, config_path=config_path)
        ConfigReader.__init__(self, config_path=config_path)
        GeometryMixin.__init__(self)
        if len(self.outlier_corrected_paths) == 0:
            raise NoFilesFoundError(msg=f'No files found in {self.outlier_corrected_dir}')
        self.session_timer = SimbaTimer(start=True)
        self.in_headers = ["Nose_x", "Nose_y", "Nose_p",
                           "Ear_left_x", "Ear_left_y", "Ear_left_p",
                           "Ear_right_x", "Ear_right_y", "Ear_right_p",
                           "Head_base_x", "Head_base_y", "Head_base_p",
                           "Tail_base_x", "Tail_base_y", "Tail_base_p"]
        self.mouse_p_headers = [x for x in self.in_headers if x[-2:] == "_p"]
        self.mouse_headers = [x for x in self.in_headers if x[-2:] != "_p"]

        
    def run(self):
        print("running")
        core_cnt, _ = find_core_cnt()
        print(f'Core count: {core_cnt}')

        for file_cnt, file_path in enumerate(self.outlier_corrected_paths):
            video_timer = SimbaTimer(start=True)
            _, video_name, _ = get_fn_ext(filepath=file_path)
            print(f'Processing {video_name}...')
            _, self.px_per_mm, fps = self.read_video_info(video_name=video_name)

            roll_windows = []
            for window in self.roll_windows_values:
                roll_windows.append(int(fps / window))

            data_df = read_df(file_path=file_path, file_type=self.file_type)
            data_df = data_df.apply(pd.to_numeric, errors='coerce').fillna(0)
            self.out_data = deepcopy(data_df)
            self.in_data = data_df
            self.in_data = self.create_shifted_df(df=self.out_data)

            save_path = os.path.join(self.features_dir, f'{video_name}.{self.file_type}')

            print("Calculating euclidean distances...")
            self.out_data["Mouse_nose_to_tail"] = self.euclidean_distance(
                self.in_data["Nose_x"].values,
                self.in_data["Tail_base_x"].values,
                self.in_data["Nose_y"].values,
                self.in_data["Tail_base_y"].values,
                self.px_per_mm,
            )
            self.out_data["Mouse_head_to_tail"] = self.euclidean_distance(
                self.in_data["Head_base_x"].values,
                self.in_data["Tail_base_x"].values,
                self.in_data["Head_base_y"].values,
                self.in_data["Tail_base_y"].values,
                self.px_per_mm,
            )

            self.out_data["Mouse_Ear_distance"] = self.euclidean_distance(
                self.in_data["Ear_left_x"].values,
                self.in_data["Ear_right_x"].values,
                self.in_data["Ear_left_y"].values,
                self.in_data["Ear_right_y"].values,
                self.px_per_mm,
            )
            self.out_data["Movement_mouse_nose"] = self.euclidean_distance(
                self.in_data["Nose_x_shifted"].values,
                self.in_data["Nose_x"].values,
                self.in_data["Nose_y_shifted"].values,
                self.in_data["Nose_y"].values,
                self.px_per_mm,
            )
            self.out_data["Movement_mouse_tail_base"] = self.euclidean_distance(
                self.in_data["Tail_base_x_shifted"].values,
                self.in_data["Tail_base_x"].values,
                self.in_data["Tail_base_y_shifted"].values,
                self.in_data["Tail_base_y"].values,
                self.px_per_mm,
            )
            self.out_data["Movement_mouse_left_ear"] = self.euclidean_distance(
                self.in_data["Ear_left_x_shifted"].values,
                self.in_data["Ear_left_x"].values,
                self.in_data["Ear_left_y_shifted"].values,
                self.in_data["Ear_left_y"].values,
                self.px_per_mm,
            )
            self.out_data["Movement_mouse_right_ear"] = self.euclidean_distance(
                self.in_data["Ear_right_x_shifted"].values,
                self.in_data["Ear_right_x"].values,
                self.in_data["Ear_right_y_shifted"].values,
                self.in_data["Ear_right_y"].values,
                self.px_per_mm,
            )
            self.out_data["Movement_mouse_head_base"] = self.euclidean_distance(
                self.in_data["Head_base_x_shifted"].values,
                self.in_data["Head_base_x"].values,
                self.in_data["Head_base_y_shifted"].values,
                self.in_data["Head_base_y"].values,
                self.px_per_mm,
            )


            print("Calculating hull variables...")
            mouse_array = self.in_data[self.mouse_headers].to_numpy()
            self.hull_dict = defaultdict(list)
            for cnt, animal_frm in enumerate(mouse_array):
                animal_frm = np.reshape(animal_frm, (-1, 2))
                animal_dists = scipy.spatial.distance.cdist(
                    animal_frm, animal_frm, metric="euclidean"
                )
                animal_dists = animal_dists[animal_dists != 0]
                self.hull_dict["M1_largest_euclidean_distance_hull"].append(
                    np.amax(animal_dists, initial=0) / self.px_per_mm
                )
                self.hull_dict["M1_smallest_euclidean_distance_hull"].append(
                    np.amin(
                        animal_dists,
                        initial=self.hull_dict["M1_largest_euclidean_distance_hull"][
                            -1
                        ],
                    )
                    / self.px_per_mm
                )
                self.hull_dict["M1_mean_euclidean_distance_hull"].append(
                    np.mean(animal_dists) / self.px_per_mm
                )
                self.hull_dict["M1_sum_euclidean_distance_hull"].append(
                    np.sum(animal_dists) / self.px_per_mm
                )
            for k, v in self.hull_dict.items():
                self.out_data[k] = v
            self.out_data["Total_movement_all_bodyparts_M1"] = (
                self.out_data["Movement_mouse_nose"]
                + self.out_data["Movement_mouse_tail_base"]
                + self.out_data["Movement_mouse_left_ear"]
                + self.out_data["Movement_mouse_right_ear"]
                + self.out_data["Movement_mouse_head_base"]
            )

            print("Calculating rolling windows: medians, medians, and sums...")

            for window in self.roll_windows_values:
                col_name = "Mouse1_mean_euclid_distances_median_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_mean_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Mouse1_mean_euclid_distances_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_mean_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Mouse1_mean_euclid_distances_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_mean_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )
            for window in self.roll_windows_values:
                col_name = "Mouse1_smallest_euclid_distances_median_{}".format(
                    str(window)
                )
                self.out_data[col_name] = (
                    self.out_data["M1_smallest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Mouse1_smallest_euclid_distances_mean_{}".format(
                    str(window)
                )
                self.out_data[col_name] = (
                    self.out_data["M1_smallest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Mouse1_smallest_euclid_distances_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_smallest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_largest_euclid_distances_median_{}".format(
                    str(window)
                )
                self.out_data[col_name] = (
                    self.out_data["M1_largest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Mouse1_largest_euclid_distances_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_largest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Mouse1_largest_euclid_distances_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["M1_largest_euclidean_distance_hull"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            for window in self.roll_windows_values:
                col_name = "Tail_base_movement_M1_median_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_tail_base"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Tail_base_movement_M1_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_tail_base"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Tail_base_movement_M1_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_tail_base"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            for window in self.roll_windows_values:
                col_name = "Head_base_movement_M1_median_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_head_base"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Head_base_movement_M1_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_head_base"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Head_base_movement_M1_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_head_base"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            for window in self.roll_windows_values:
                col_name = "Nose_movement_M1_median_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_nose"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Nose_movement_M1_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_nose"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Nose_movement_M1_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Movement_mouse_nose"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            for window in self.roll_windows_values:
                col_name = "Total_movement_M1_median_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Total_movement_all_bodyparts_M1"]
                    .rolling(int(window), min_periods=1)
                    .median()
                )
                col_name = "Total_movement_M1_mean_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Total_movement_all_bodyparts_M1"]
                    .rolling(int(window), min_periods=1)
                    .mean()
                )
                col_name = "Total_movement_M1_sum_{}".format(str(window))
                self.out_data[col_name] = (
                    self.out_data["Total_movement_all_bodyparts_M1"]
                    .rolling(int(window), min_periods=1)
                    .sum()
                )

            print("Calculating deviations...")
            self.out_data["Total_movement_all_bodyparts_deviation"] = (
                self.out_data["Total_movement_all_bodyparts_M1"].mean()
                - self.out_data["Total_movement_all_bodyparts_M1"]
            )
            self.out_data["M1_smallest_euclid_distances_hull_deviation"] = (
                self.out_data["M1_smallest_euclidean_distance_hull"].mean()
                - self.out_data["M1_smallest_euclidean_distance_hull"]
            )
            self.out_data["M1_largest_euclid_distances_hull_deviation"] = (
                self.out_data["M1_largest_euclidean_distance_hull"].mean()
                - self.out_data["M1_largest_euclidean_distance_hull"]
            )
            self.out_data["M1_mean_euclid_distances_hull_deviation"] = (
                self.out_data["M1_mean_euclidean_distance_hull"].mean()
                - self.out_data["M1_mean_euclidean_distance_hull"]
            )

            for window in self.roll_windows_values:
                col_name = "Mouse1_smallest_euclid_distances_mean_{}".format(
                    str(window)
                )
                deviation_col_name = col_name + "_deviation"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_largest_euclid_distances_mean_{}".format(str(window))
                deviation_col_name = col_name + "_deviation"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_mean_euclid_distances_mean_{}".format(str(window))
                deviation_col_name = col_name + "_deviation"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            '''
            Old script:
            self.out_data["Movement_mouse_nose_percentile_rank"] = self.out_data[
                "Movement_mouse_nose"
            ].rank(pct=True)
            '''

            print("Calculating percentile ranks...")
            # List of body parts to process
            body_parts = ["mouse_nose", "mouse_head_base"]
            # Loop through each body part and calculate percentile rank
            for body_part in body_parts:
                col_name = f"Movement_{body_part}"  # Dynamically create the column name
                percentile_col_name = f"Movement_{body_part}_percentile_rank"  # Output column name
                self.out_data[percentile_col_name] = self.out_data[col_name].rank(pct=True)

            for window in self.roll_windows_values:
                col_name = "Total_movement_M1_mean_{}".format(str(window))
                deviation_col_name = col_name + "_percentile_rank"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_mean_euclid_distances_mean_{}".format(str(window))
                deviation_col_name = col_name + "_percentile_rank"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_smallest_euclid_distances_mean_{}".format(
                    str(window)
                )
                deviation_col_name = col_name + "_percentile_rank"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            for window in self.roll_windows_values:
                col_name = "Mouse1_largest_euclid_distances_mean_{}".format(str(window))
                deviation_col_name = col_name + "_percentile_rank"
                self.out_data[deviation_col_name] = (
                    self.out_data[col_name].mean() - self.out_data[col_name]
                )

            print("Calculating path tortuosities...")
            as_strided = np.lib.stride_tricks.as_strided
            win_size = 3
            centroid_lst_mouse_x = as_strided(
                self.out_data["Nose_x"],
                (len(self.out_data) - (win_size - 1), win_size),
                (self.out_data["Nose_x"].values.strides * 2),
            )
            centroid_lst_mouse_y = as_strided(
                self.out_data["Nose_y"],
                (len(self.out_data) - (win_size - 1), win_size),
                (self.out_data["Nose_y"].values.strides * 2),
            )

            for window in self.roll_windows_values:
                start, end = 0, 0 + int(window)
                tortuosities_results = defaultdict(list)
                for frame in range(len(self.out_data)):
                    tortuosities_dict = defaultdict(list)
                    c_centroid_lst_mouse_x, c_centroid_lst_mouse_y = (
                        centroid_lst_mouse_x[start:end],
                        centroid_lst_mouse_y[start:end],
                    )
                    for frame_in_window in range(len(c_centroid_lst_mouse_x)):
                        move_angle_mouse_ = self.angle3pt(
                            c_centroid_lst_mouse_x[frame_in_window][0],
                            c_centroid_lst_mouse_y[frame_in_window][0],
                            c_centroid_lst_mouse_x[frame_in_window][1],
                            c_centroid_lst_mouse_y[frame_in_window][1],
                            c_centroid_lst_mouse_x[frame_in_window][2],
                            c_centroid_lst_mouse_y[frame_in_window][2],
                        )
                        tortuosities_dict["Animal_1"].append(move_angle_mouse_)
                    tortuosities_results["Animal_1"].append(
                        sum(tortuosities_dict["Animal_1"]) / (2 * math.pi)
                    )
                    start += 1
                    end += 1
                col_name = "Tortuosity_Mouse1_{}".format(str(window))
                self.out_data[col_name] = tortuosities_results["Animal_1"]

            print("Calculating pose probability scores...")
            self.out_data["Sum_probabilities"] = (
                self.out_data["Ear_left_p"]
                + self.out_data["Ear_right_p"]
                + self.out_data["Nose_p"]
                + self.out_data["Tail_base_p"]
                + self.out_data["Head_base_p"]
            )
            self.out_data["Sum_probabilities_deviation"] = (
                self.out_data["Sum_probabilities"].mean()
                - self.out_data["Sum_probabilities"]
            )
            self.out_data["Sum_probabilities_deviation_percentile_rank"] = (
                self.out_data["Sum_probabilities_deviation"].rank(pct=True)
            )
            self.out_data["Sum_probabilities_percentile_rank"] = self.out_data[
                "Sum_probabilities_deviation_percentile_rank"
            ].rank(pct=True)
            results = pd.DataFrame(
                self.count_values_in_range(
                    data=self.out_data.filter(self.mouse_p_headers).values,
                    ranges=np.array([[0.0, 0.1], [0.0, 0.5], [0.0, 0.75]]),
                ),
                columns=[
                    "Low_prob_detections_0.1",
                    "Low_prob_detections_0.5",
                    "Low_prob_detections_0.75",
                ],
            )


            print(f"We now have {len(self.out_data.columns)} columns")
            self.save(data=self.out_data, save_path=save_path)
            video_timer.stop_timer()
            stdout_success(msg=f'{video_name} complete!', elapsed_time=video_timer.elapsed_time_str)

    def save(self,
             data: pd.DataFrame,
             save_path: str):

        write_df(df=data, file_type=self.file_type, save_path=save_path)