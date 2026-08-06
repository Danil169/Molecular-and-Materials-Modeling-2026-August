from pathlib import Path

def list_files_by_size(directory_path, descending=True):
    # Convert string path to a Path object
    path = Path(directory_path)
    
    # Check if the path exists and is a directory
    if not path.is_dir():
        print(f"Error: {directory_path} is not a valid directory.")
        return

    # Filter out directories, keeping only files
    files = [item for item in path.iterdir() if item.is_file()]
    
    # Sort files by their size attribute (stat().st_size)
    files.sort(key=lambda file: file.stat().st_size, reverse=descending)
    
    # Print the sorted results
    print(f"{'File Name':<40} | {'Size (Bytes)':<12}")
    print("-" * 55)
    for file in files:
        print(f"{file.name:<40} | {file.stat().st_size:<12,}")

# Example Usage: Replace '.' with your target folder path (e.g., "C:/Users/Name/Documents")
list_files_by_size('.')

