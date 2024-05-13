import os
import re
import subprocess

def find_python_imports(root_dir):
    import_statements = set()
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(subdir, file)
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        contents = f.read()
                except UnicodeDecodeError:
                    # Try reading with 'latin-1' encoding if UTF-8 fails
                    with open(file_path, 'r', encoding='latin-1') as f:
                        contents = f.read()
                imports = re.findall(r'^import (\S+)|^from (\S+) import', contents, re.MULTILINE)
                for imp in imports:
                    # Add the first non-empty matching group
                    import_statements.add(next(filter(None, imp)))
    return import_statements

def execute_pyinstaller(main_script, hidden_imports):
    cmd = ["pyinstaller", "--onedir"]
    cmd.extend(f"--hidden-import={imp}" for imp in hidden_imports)
    cmd.append(main_script)
    subprocess.run(cmd)

# Replace 'path/to/your/project' with the actual path to your project
project_path = './'
# Replace 'your_script.py' with the name of your main Python script
main_script = 'satute_cli.py'

imports = find_python_imports(project_path)
execute_pyinstaller(main_script, imports)

print("PyInstaller command executed with --onedir. Check the 'dist' directory for your application.")
