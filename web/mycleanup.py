import os
import time

directory = 'userFiles/'
Min=(100,86400)
Moderate=(1000,3600)
Greedy=(5000,600)

current_time = time.time()
files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
file_count = len(files)

try:
    files.remove("index.html")
except ValueError:
    pass  # 'index.html' was not in the list

# Determine the age limit based on the file count
if file_count > Min[0]:
    if file_count>Greedy[0]:
        max_age=Greedy[1]
    elif file_count>Moderate[0]:
        max_age=Moderate[1]
    else:
        max_age=Min[1]

    # Delete files older than the determined max_age
    for file in files:
        file_path = os.path.join(directory, file)
        file_creation_time = os.path.getmtime(file_path)
        if current_time - file_creation_time > max_age:
            try: 
                os.remove(file_path)
            except:
                print(file_path)


