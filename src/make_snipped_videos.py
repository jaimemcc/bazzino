#!/usr/bin/env python
"""
Extract and save video snippets around TTL events.

This script processes behavioral video data aligned with TTL (transistor-transistor logic) 
events from TDT (Tucker-Davis Technologies) recordings. It extracts video snippets around 
each event and saves them with optional text overlays.

Usage:
    python make_snipped_videos.py --stub PB71-221123-113609 --ttls-csv ../data/ttls.csv -i D://TestData//bazzino//from_paula -o ../results
"""

import argparse
import cv2
from pathlib import Path
import pandas as pd


def load_ttls_from_csv(csv_path, stub):
    """
    Load TTL times from a CSV file for a specific recording.
    
    Parameters:
    -----------
    csv_path : Path
        Path to CSV file containing TTL data with stub as a column
    stub : str
        Recording identifier (column name in CSV)
        
    Returns:
    --------
    np.ndarray
        Array of TTL times in seconds
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"TTL CSV file not found: {csv_path}")
    
    df = pd.read_csv(csv_path, index_col=0)
    
    if stub not in df.columns:
        raise ValueError(f"Recording '{stub}' not found in TTL CSV. Available recordings: {list(df.columns)}")
    
    # Remove NaN values and return as numpy array
    ttls = df[stub].dropna().values
    return ttls


def get_videopath(stub):
    """
    Construct input video filename from recording stub.
    
    Parameters:
    -----------
    stub : str
        Recording identifier (e.g., 'PB71-221123-113609')
        
    Returns:
    --------
    str
        Filename of the video file
    """
    date = stub.split("-")[1]
    return f"PB_NAapp-{date}_{stub}_Cam1.avi"


def save_video_snippet(input_path, output_path, start_frame, end_frame, 
                       stub=None, trial_num=None, text=None, text_start_frame=50):
    """
    Extract a portion of a video file and save it.
    
    Parameters:
    -----------
    input_path : Path or str
        Path to input video file
    output_path : Path or str
        Root output directory where organized output will be saved
    start_frame : int
        Starting frame number (0-indexed)
    end_frame : int
        Ending frame number (exclusive)
    stub : str, optional
        Recording identifier - used to organize output in subdirectories
    trial_num : int, optional
        Trial number to include in output filename
    text : str, optional
        Text to overlay on the video
    text_start_frame : int, default=50
        Frame number (relative to start_frame) to start displaying text
    """
    # Convert Path objects
    input_path = Path(input_path)
    output_path = Path(output_path)
    
    # Detect input video format from extension
    input_ext = input_path.suffix.lower()
    
    # Check if input video is filtered
    is_filtered = "filtered" in input_path.name.lower()
    filtered_suffix = "_filtered" if is_filtered else ""
    
    # Create output directory if using stub
    if stub:
        output_folder = output_path / stub
        output_folder.mkdir(parents=True, exist_ok=True)
        if trial_num is not None:
            output_file = output_folder / f"{stub}_trial_{trial_num:02d}_{start_frame}{filtered_suffix}{input_ext}"
        else:
            output_file = output_folder / f"{stub}_{start_frame}{filtered_suffix}{input_ext}"
    else:
        output_folder = output_path
        output_folder.mkdir(parents=True, exist_ok=True)
        output_file = output_folder / f"snippet_{start_frame}_{end_frame}{filtered_suffix}{input_ext}"
    
    output_path_str = str(output_file)
    print(f"Saving video to: {output_path_str}")

    # Open the input video
    cap = cv2.VideoCapture(str(input_path))
    
    if not cap.isOpened():
        print(f"Error: Could not open input video: {input_path}")
        return
    
    # Get video properties
    fps = cap.get(cv2.CAP_PROP_FPS)
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    
    print(f"Input video properties: {width}x{height} @ {fps} fps")
    
    # Select appropriate codec based on output format
    out = None
    if input_ext == '.mp4':
        # Try H.264 codec for MP4
        codec_options = ['mp4v', 'avc1', 'H264']
    else:
        # AVI format codecs
        codec_options = ['XVID', 'MJPG']
    
    for codec in codec_options:
        fourcc = cv2.VideoWriter_fourcc(*codec)
        out = cv2.VideoWriter(output_path_str, fourcc, fps, (width, height))
        if out.isOpened():
            print(f"Using codec: {codec}")
            break
        else:
            print(f"Warning: Could not create VideoWriter with {codec}. Trying next codec...")
    
    if out is None or not out.isOpened():
        print(f"Error: Could not create VideoWriter with any available codec")
        cap.release()
        return
    
    # Set starting frame
    cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)
    
    # Extract and write frames
    frame_count = start_frame
    snippet_frame = 0
    while frame_count < end_frame:
        ret, frame = cap.read()
        if not ret:
            break
        
        # Add text if specified and frame count exceeds text_start_frame
        if text is not None and snippet_frame >= text_start_frame:
            cv2.putText(frame, text, (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 
                        1, (0, 255, 0), 2)
        
        out.write(frame)
        frame_count += 1
        snippet_frame += 1
    
    # Release everything
    cap.release()
    out.release()
    
    num_frames = end_frame - start_frame
    print(f"Video snippet saved: {num_frames} frames")


def extract_videos(stub, input_folder, output_folder, ttls_csv, pre_seconds=5, post_seconds=15):
    """
    Extract video snippets for all TTL events in a recording.
    
    Parameters:
    -----------
    stub : str
        Recording identifier
    input_folder : Path
        Path to a video file or root folder containing video files
    output_folder : Path
        Output folder for video snippets
    ttls_csv : Path
        Path to CSV file containing TTL data
    pre_seconds : float, default=5
        Seconds of video to extract before each TTL event
    post_seconds : float, default=15
        Seconds of video to extract after each TTL event
    """
    ttls = load_ttls_from_csv(ttls_csv, stub)
    
    if len(ttls) == 0:
        print(f"No TTLs found for stub {stub}")
        return
    
    # Check if input_folder is actually a video file
    input_path = Path(input_folder)
    video_extensions = ('.avi', '.mp4', '.mov', '.mkv', '.flv', '.wmv', '.webm', '.m4v')
    
    if input_path.is_file() and input_path.suffix.lower() in video_extensions:
        video_path = input_path
    else:
        # Treat as a folder and construct the video path from stub
        video_path = input_path / get_videopath(stub)
        
        if not video_path.exists():
            video_path = input_path / stub / get_videopath(stub)
        
        if not video_path.exists():
            raise FileNotFoundError(f"Video file not found for stub {stub} at expected locations: {video_path}")
    
    # Get video properties to determine frame rate
    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        print(f"Error: Could not open video: {video_path}")
        return
    
    fps = cap.get(cv2.CAP_PROP_FPS)
    cap.release()
    
    print(f"Processing {len(ttls)} TTL events from {stub}")
    
    # Extract snippet for each TTL
    for i, ttl in enumerate(ttls):
        start = int((ttl - pre_seconds) * fps)
        end = int((ttl + post_seconds) * fps)
        
        # Ensure frames are within valid range
        start = max(0, start)
        
        print(f"  TTL {i+1}/{len(ttls)}: frames {start}-{end}")
        save_video_snippet(video_path, output_folder, 
                          start_frame=start, end_frame=end, 
                          stub=stub, trial_num=i, text="Infusion")


def main():
    parser = argparse.ArgumentParser(
        description="Extract and save video snippets around TTL events from CSV data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python make_snipped_videos.py --stub PB71-221123-113609 --ttls-csv ../data/ttls.csv \\
    -i D://TestData//bazzino//from_paula -o ../results
        """)
    
    parser.add_argument("--stub", type=str, required=True,
                       help="Recording identifier")
    parser.add_argument("--ttls-csv", type=Path, required=True,
                       help="Path to CSV file containing TTL data")
    parser.add_argument("-i", "--input_folder", type=Path, required=True, dest="input_folder",
                       help="Path to a video file or root folder containing video files")
    parser.add_argument("-o", "--output_folder", type=Path, required=True, dest="output_folder",
                       help="Output folder for video snippets")
    parser.add_argument("--pre", type=float, default=5,
                       help="Seconds of video to extract before TTL event (default: 5)")
    parser.add_argument("--post", type=float, default=15,
                       help="Seconds of video to extract after TTL event (default: 15)")
    args = parser.parse_args()
    
    try:
        extract_videos(args.stub, args.input_folder, args.output_folder, 
                      ttls_csv=args.ttls_csv, pre_seconds=args.pre, post_seconds=args.post)
        print("Video extraction complete")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        parser.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        parser.exit(1)


if __name__ == "__main__":
    main()

