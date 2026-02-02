from pathlib import Path

def find_folders_without_avis(root_folder):
    """
    Find all folders within root_folder that do NOT contain any AVI files.
    
    Parameters:
    - root_folder: path to the root folder to search
    
    Returns:
    - List of folder paths that don't contain AVI files
    """
    root = Path(root_folder)
    
    if not root.exists():
        print(f"Error: Folder does not exist: {root_folder}")
        return []
    
    folders_without_avis = []
    
    # Get all subdirectories
    for folder in root.rglob('*'):
        if folder.is_dir():
            # Check if this folder contains any .avi files
            avi_files = list(folder.glob('*.avi'))
            
            if not avi_files:
                folders_without_avis.append(folder)
    
    return folders_without_avis

if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        root_folder = sys.argv[1]
    else:
        root_folder = input("Enter the root folder path: ")
    
    result = find_folders_without_avis(root_folder)
    
    print(f"\nFolders without AVI files ({len(result)} total):\n")
    for folder in sorted(result):
        print(folder)
    
    # Save to a text file
    output_file = Path(root_folder) / "folders_without_avis.txt"
    with open(output_file, 'w') as f:
        for folder in sorted(result):
            f.write(f"{folder}\n")
    
    print(f"\nResults saved to: {output_file}")
