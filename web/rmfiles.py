import os
import sys

def main(accode):
    directory = 'userFiles/'
    files = os.listdir(directory)
    for file in files:
        # Check if 'accode' is in the filename and it does not start with "Input"
        if accode in file and not file.startswith("input"):
            try:
                os.remove(os.path.join(directory, file))
                print(f"Deleted: {file}")
            except Exception as e:
                print(f"Failed to delete {file}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python removefiles.py <accode>")
        sys.exit(1)
    accode = sys.argv[1]
    main(accode)
