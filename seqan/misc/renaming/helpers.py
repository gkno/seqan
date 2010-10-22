import os
import os.path

def test_file_type(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    return ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++"]

def all_files(project_path):
    for root, dirs, files in os.walk(project_path):
        if 'CVS' in dirs:
            dirs.remove('CVS')
        if '.svn' in dirs:
            dirs.remove('.svn')
        for file in files:
            if file.startswith('.'):
                continue # Skip hidden files.
            if file.find('_generated_forwards') != -1:
                continue # Skip generated forwards.
            path = os.path.join(root, file)
            if test_file_type(path):
                yield path
