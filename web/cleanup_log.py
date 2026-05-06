import os
import time

directory = 'log/'
maxn=250

current_time = time.time()
files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
file_count = len(files)

try:
    files.remove("index.html")
except ValueError:
    pass  # 'index.html' was not in the list

# Remove oldest files if file size is greater than maxn
if file_count - 1 > maxn:
    full_paths = [os.path.join(directory, f) for f in files]
    full_paths.sort(key=os.path.getmtime)  # oldest first
    files_to_delete = full_paths[:file_count - maxn]
    for f in files_to_delete:
        os.remove(f)




