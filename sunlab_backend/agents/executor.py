import subprocess

class Executor:
    def __init__(self, llm, data_path, memory):
        self.llm = llm
        self.data_path = data_path
        self.memory = memory

    def run_python_code(self, code: str):
        with open("temp_code.py", "w") as f:
            f.write(code)
        try:
            result = subprocess.check_output(["python", "temp_code.py"], stderr=subprocess.STDOUT)
            return result.decode()
        except subprocess.CalledProcessError as e:
            return e.output.decode()
